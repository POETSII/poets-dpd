#ifndef naive_dpd_engine_half_step_hpp
#define naive_dpd_engine_half_step_hpp

#include "dpd/core/dpd_engine.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"

#include <cassert>
#include <cmath>
#include <array>

/*
    This method must have lambda=0.5

    The process to get from t to t+dt is:
    - UpdatePos : calculate position at time t+dt using values only from time t
        // Pre: x(t),v(t),f(t) all live
        x(t+dt)   = x(t) + v(t)*dt + f(t)*(dt^2/2)
        v(t+dt/2) = v(t) + f(t)*(dt/2),  This assumes lambda=0.5
        // Post: x(t+dt),v(t+dt) live,   x(t) is dead, f(t) is dead

    - UpdateForce : calculate forces using positions at time t+dt and velocities at t+dt/2
        // Pre: x(t+dt),v(t+dt/2)
        f(t+dt) = (function of x(t+dt),v(t+dt/2))
        // Post: f(t+dt),x(t+dt),f(t),v(t) all live

    - UpdateMom : calculate velocity at time t+dt
        v(t+dt) = v(t+dt/2) + dt*f(t+dt)/2
        // Post: f(t+dt),v(t+dt),x(t+dt) all live, v(t+dt/2) dead

    To keep state clean, the approach taken here is not to pre-calculate
    v'(t), and just to recalculate on the fly. This engine is already
    pretty inefficient anyway.
*/
class NaiveDPDEngineHalfStep
    : public DPDEngine
{
    friend class NaiveDPDEngineHalfStepTBB;
    friend class NaiveDPDEngineHalfStepTBBV2;
    
public:
    virtual void Attach(WorldState *state)
    {
        m_state=state;
    }

    virtual void Run(unsigned nSteps) override
    {
        check_constraints_and_setup();

        for(unsigned i=0; i<nSteps; i++){
            step();
        }
    }

private:
    WorldState *m_state;

    unsigned m_numBeadTypes;
    vec3r_t m_origins; // origin in each dimension. Must be an integer, though encoded as real
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    std::vector<std::vector<Bead*>> m_cells;

    double m_inv_root_dt;

    uint64_t m_t_hash;

    void check_constraints_and_setup()
    {
        auto require=[](bool cond, const char *msg)
        {
            if(!cond){
                throw std::string(msg);
            }
        };

        require( m_state, "No state" );
        for(unsigned i=0; i<3; i++){
            require( round(m_state->box[i]) == m_state->box[i], "box must be integer aligned");
            require( m_state->box[i] >= 2, "Distance in each direction must be at least 2");
            m_lengths[i]=m_state->box[i];
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_inv_root_dt=recip_pow_half(m_state->dt);
    }

    unsigned calc_num_cells() const
    {
        unsigned n=1;
        for(unsigned i=0;i<3;i++){
            n *= m_lengths[i];
        }
        return n;
    }

    vec3i_t world_pos_to_cell_pos(const vec3r_t &pos) const
    {
        vec3i_t res;
        for(unsigned i=0; i<3; i++){
            res[i] = floor(pos[i] - m_origins[i]);
            assert(0<=res[i] && res[i] < m_lengths[i]);
        }
        return res;
    }

    unsigned cell_pos_to_index(vec3i_t pos) const
    {
        unsigned index=0;
        unsigned last_dim=0;
        for(unsigned i=0; i<3; i++){
            assert(0<= pos[i] && pos[i] < m_lengths[i] );
            index = index * last_dim + pos[i];
            last_dim = m_lengths[i];
        }
        return index;
    }

    unsigned world_pos_to_cell_index(const vec3r_t &pos) const
    { return cell_pos_to_index(world_pos_to_cell_pos(pos)); }

    // This produces a pair (cell_loc,pos_adj)
    // - cell_loc : the cell position of the cell in direction dir
    // - delta : the delta to add to beads from the cell in direction dir to get correct relative position
    std::pair<vec3i_t,vec3r_t> make_relative_cell_pos(const vec3i_t &base, const vec3i_t &dir) const
    {
        vec3i_t pos;
        vec3r_t delta;
        for(unsigned i=0; i<3; i++){
            assert(-1 <= dir[i] && dir[i] <= +1);
            int raw = base[i] + dir[i];
            int base_is_left = raw < 0;
            int base_is_right = raw == m_lengths[i];

            pos[i] = raw + (base_is_left - base_is_right) * m_lengths[i];

            delta[i] = (base_is_right - base_is_left) * m_lengths[i];
        }
        return {pos,delta};
    }

    vec3r_t calc_distance_from_to(const vec3r_t &base, const vec3r_t &other) const
    {
        vec3r_t res;
        for(unsigned i=0; i<3; i++){
            double len=m_lengths[i];
            double dx=other[i]-base[i];
            int wrap_neg = dx > len*0.5;
            int wrap_pos = dx < -len*0.5;
            res[i] = dx + (wrap_pos - wrap_neg) * len;
        }
        return res;
    }

    virtual void step()
    {
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        double dt=m_state->dt;

        // Clear all cell information
        m_cells.resize(calc_num_cells());
        for(auto &c : m_cells){
            c.clear();
        }

        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &b : m_state->beads){
            dpd_maths_core_half_step::update_pos(dt, m_state->box, b);
            m_cells.at( world_pos_to_cell_index(b.x) ).push_back(&b);
        }

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        vec3i_t pos, dir;
        for(pos[0]=0; pos[0]<m_lengths[0]; pos[0]++){
            for(pos[1]=0; pos[1]<m_lengths[1]; pos[1]++){
                for(pos[2]=0; pos[2]<m_lengths[2]; pos[2]++){
                    auto &home=m_cells[cell_pos_to_index(pos)];
                    for(dir[0]=-1; dir[0]<=+1; dir[0]++){
                        for(dir[1]=-1; dir[1]<=+1; dir[1]++){
                            for(dir[2]=-1; dir[2]<=+1; dir[2]++){
                                auto [other_pos,other_delta] = make_relative_cell_pos(pos, dir);
                                const auto &other=m_cells[cell_pos_to_index(other_pos)];
                                update_cell_forces(home, other, other_delta);
                            }
                        }
                    }
                }
            }
        }

        // Update all angle bonds
        for(const auto &p : m_state->polymers){
            const auto &pt = m_state->polymer_types.at(p.polymer_type);
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(p, pt, bond_pair);
            }
        }

        // Final momentum update
        for(auto &b : m_state->beads){
            dpd_maths_core_half_step::update_mom(dt, b);
        }

        m_state->t += 1;
    }

    bool is_bonded(
        const Bead &home, const Bead &other,
        double &kappa,
        double &r0
    ) const
    {
        kappa=0;
        r0=0;
        if(home.polymer_id!=other.polymer_id){
            return false; // Rejects the vast majority
        }
        // This is very slow. In practise use a better method
        for(const Bond &b : m_state->polymer_types.at(home.polymer_type).bonds){
            if( ((b.bead_offset_head==home.polymer_offset) && (b.bead_offset_tail==other.polymer_offset))
                ||
                ((b.bead_offset_head==other.polymer_offset) && (b.bead_offset_tail==home.polymer_offset))
            ){
                kappa=b.kappa;
                r0=b.r0;
                return true;
            }
        }
        return false;
    }

    template<class TBeadPtrVec>
    void update_cell_forces(TBeadPtrVec &home, const TBeadPtrVec &other, const vec3r_t &other_delta) const
    {
        for(Bead *hb : home){
            for(const Bead *ob : other)
            {
                if(hb==ob){
                    continue;
                }

                vec3r_t dx = vec3r_t(hb->x) - vec3r_t(ob->x) - other_delta; 
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<0.00001){
                    continue;
                }

                double dr=pow_half(dr2);

                double kappa,r0;
                is_bonded(*hb, *ob, kappa, r0);

                vec3r_t f;
                
                dpd_maths_core_half_step::calc_force(
                    (recip_pow_half(m_state->dt)),
                    [&](unsigned a, unsigned b){ return m_state->interactions[a*m_numBeadTypes+b].conservative; },
                    [&](unsigned a, unsigned b){ return pow_half(m_state->interactions[a*m_numBeadTypes+b].dissipative); },
                    m_t_hash,
                    dx, dr,
                    kappa, r0,
                    *hb,
                    *ob,
                    f
                );
                hb->f += f;
            }
        }
    }

    void update_angle_bond(const Polymer &p, const PolymerType &pt, const BondPair &bp) const
    {
        unsigned head_bond_index=bp.bond_offset_head;
        unsigned tail_bond_index=bp.bond_offset_tail;
        const Bond &head_bond=pt.bonds.at(head_bond_index);
        const Bond &tail_bond=pt.bonds.at(tail_bond_index);

        unsigned head_bead_index=head_bond.bead_offset_head;
        unsigned middle_bead_index=head_bond.bead_offset_tail;
        assert(middle_bead_index==tail_bond.bead_offset_head);
        unsigned tail_bead_index=tail_bond.bead_offset_tail;

        Bead &head_bead=m_state->beads.at(p.bead_ids.at(head_bead_index));
        Bead &middle_bead=m_state->beads.at(p.bead_ids.at(middle_bead_index));
        Bead &tail_bead=m_state->beads.at(p.bead_ids.at(tail_bead_index));

        auto first=calc_distance_from_to(head_bead.x, middle_bead.x);
        auto second=calc_distance_from_to(middle_bead.x, tail_bead.x);

        double FirstLength   = first.l2_norm();
        double SecondLength  = second.l2_norm();

        vec3r_t headForce, middleForce, tailForce;

        dpd_maths_core_half_step::calc_angle_force(
            bp.kappa, cos(bp.theta0), sin(bp.theta0), 
            first, FirstLength,
            second, SecondLength,
            headForce, middleForce, tailForce
        );
        //std::cerr<<"Dur: dx01="<<first<<", dx12="<<second<<", headForce="<<headForce<<"\n";

        head_bead.f += headForce;
        tail_bead.f += tailForce;
        middle_bead.f += middleForce;
    }
};

#endif