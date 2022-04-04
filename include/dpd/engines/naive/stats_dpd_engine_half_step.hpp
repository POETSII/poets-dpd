#ifndef stats_dpd_engine_half_step_hpp
#define stats_dpd_engine_half_step_hpp

#include "dpd/core/dpd_engine.hpp"
#include "dpd/engines/naive/naive_dpd_engine.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"

#include <cassert>
#include <cmath>
#include <array>
#include <unordered_set>

#include <immintrin.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/min.hpp>

#include <cmath>

namespace dpd_stats
{
    using namespace boost::accumulators;

    struct count_stats
    {
        std::vector<uint64_t> counts;

        void add(unsigned i)
        {
            if(i>=counts.size()){
                counts.resize(i+1, 0);
            }
            counts[i] += 1;
        }

        void add(const count_stats &s)
        {
            if(s.counts.size() > counts.size()){
                counts.resize(s.counts.size());
            }
            for(unsigned i=0; i<s.counts.size(); i++){
                counts[i] += s.counts[i];
            }
        }
    };

    struct float_stats{
        double pos_smallest, pos_biggest;
        double neg_smallest, neg_biggest;
        

        void add(double x)
        {
            acc(x);
        }
    };
};

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
class StatsDPDEngineHalfStep
    : public DPDEngine
{
    
public:
    enum InteractionDir
    {
        Centre,
        Face,
        Edge,
        Corner
    };

    static InteractionDir calc_interaction_dir(int dx, int dy, int dz)
    {
        switch(std::abs(dx)+std::abs(dy)+std::abs(dz)){
        case 0: return Centre;
        case 1: return Face;
        case 2: return Edge;
        case 3: return Corner;
        default: assert(0); return Centre;
        }
    }

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
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    vec3i_t m_dims;// Same as m_leNgths
    std::vector<std::vector<Bead*>> m_cells;

    double m_scaled_inv_root_dt;
    uint64_t m_t_hash;

    struct interaction_stats
    {
        dpd_stats::float_stats r;
        dpd_stats::count_stats hit;
    };

    struct round_stats
    {
        unsigned num_rounds;
        dpd_stats::count_stats migrations_total;

        dpd_stats::count_stats beads_per_cell;
        dpd_stats::count_stats migrations_in_per_cell;
        dpd_stats::count_stats migrations_out_per_cell;

        interaction_stats interactions_all;
        std::array<interaction_stats,4> interactions_by_dir;
    };

    round_stats m_curr_stats;


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
            m_dims[i]=(int)m_state->box[i];
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_scaled_inv_root_dt=pow_half(24*dpd_maths_core_half_step::kT / m_state->dt);
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
            res[i] = (int)floor(pos[i]);
            assert(0<=res[i] && res[i] < m_dims[i]);
        }
        return res;
    }

    unsigned cell_pos_to_index(vec3i_t pos) const
    {
        return pos[0]*m_dims[1]*m_dims[2] + pos[1] * m_dims[2] + pos[2];
    }

    vec3i_t index_to_cell_pos(unsigned index) const
    {
        return { (int)(index/(m_dims[1]*m_dims[2])) , (int)((index/m_dims[2])%m_dims[0]) , (int)(index%m_dims[2]) };
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

            pos[i] = raw + (base_is_left - base_is_right) * m_dims[i];

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
        step_impl();
    }

    void step_impl()
    {
        m_curr_stats.num_rounds++;

        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        double dt=m_state->dt;

        // Clear all cell information
        m_cells.resize(calc_num_cells());
        for(auto &c : m_cells){
            c.clear();
        }

        std::vector<unsigned> migrate_in(m_cells.size(), 0);
        std::vector<unsigned> migrate_out(m_cells.size(), 0);

        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &b : m_state->beads){
            vec3i_t ppos=vec3_floor(b.x);
            dpd_maths_core_half_step::update_pos(dt, m_state->box, b);
            vec3i_t npos=vec3_floor(b.x);

            if(ppos!=npos){
                migrate_out[cell_pos_to_index(ppos)]+=1;
                migrate_in[cell_pos_to_index(npos)]+=1;
            }
            
            m_cells.at( world_pos_to_cell_index(b.x) ).push_back(&b);
        }

        for(unsigned i=0; i<m_cells.size(); i++){
            m_curr_stats.migrations_in_per_cell.add(migrate_in[i]);
            m_curr_stats.migrations_out_per_cell.add(migrate_out[i]);
            m_curr_stats.beads_per_cell.add(m_cells[i].size());
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
                                update_cell_forces<false>(home, other, other_delta, calc_interaction_dir(dir[0], dir[1], dir[2]));
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

    template<bool EnableLogging,class TBeadPtrVec>
    void  __attribute__((noinline))  update_cell_forces(TBeadPtrVec &home, const TBeadPtrVec &other, const vec3r_t &other_delta, InteractionDir dir)
    {        
        for(const Bead *ob : other)
        {
            vec3r_t ob_x = vec3r_t(ob->x) + other_delta;
            for(Bead *hb : home){
                vec3r_t dx =  hb->x - ob_x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                double dr=pow_half(dr2);

                bool hit = dr2 < 1 && dr2>=MIN_DISTANCE_CUTOFF_SQR;

                m_curr_stats.interactions_all.r.add(dr);
                m_curr_stats.interactions_all.hit.add(hit);
                m_curr_stats.interactions_by_dir[dir].r.add(dr);
                m_curr_stats.interactions_by_dir[dir].hit.add(hit);
                
                if(!hit){
                    continue;
                }

                assert(hb!=ob);

                double kappa,r0;
                is_bonded(*hb, *ob, kappa, r0);

                vec3r_t f;
                dpd_maths_core_half_step::calc_force<EnableLogging>(
                    m_scaled_inv_root_dt,
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

        head_bead.f += headForce;
        tail_bead.f += tailForce;
        middle_bead.f += middleForce;
    }
};

#endif
