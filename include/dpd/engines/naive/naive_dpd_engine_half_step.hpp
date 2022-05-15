#ifndef naive_dpd_engine_half_step_hpp
#define naive_dpd_engine_half_step_hpp

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
    friend class TolerantDPDEngineHalfStepTBB;
    
    static double now()
    {
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return ts.tv_sec+1e-9*ts.tv_nsec;
    };
public:
    virtual void Attach(WorldState *state)
    {
        m_state=state;
        m_timings=DPDEngine::timings_t();
    }

    virtual void Run(unsigned nSteps) override
    {
        double start=now();

        check_constraints_and_setup();

        double mid=now();
        m_timings.configure += mid-start;

        for(unsigned i=0; i<nSteps; i++){
            step();
        }

        double end=now();
        m_timings.execute_to_first_bead += end-mid;
        m_timings.execute_to_last_bead += end-mid; 
    }

    void SetFixedWater(bool fixed_water=true)
    {
        m_fixed_water=fixed_water;
    }

    // This is a bit arbitrary, but this is not supposed
    // to be an engine that operates in weird stressed situations.
    double m_max_bond_length=2.0;

    virtual double GetMaxBondLength() const
    {
        return m_max_bond_length;
    }

    virtual bool GetTimings(DPDEngine::timings_t &timings)
    {
        timings = m_timings;
        return true;
    }
private:
    WorldState *m_state;

    timings_t m_timings;

    unsigned m_numBeadTypes;
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    vec3i_t m_dims;// Same as m_leNgths
    std::vector<std::vector<Bead*>> m_cells;

    double m_scaled_inv_root_dt;

    uint64_t m_t_hash;
    bool m_fixed_water=false;

    std::unordered_set<std::pair<const Bead*,const Bead*>, NaiveDPDEngine<false>::bp_hash> m_seen_pairs;

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
        m_seen_pairs.clear();

        if(ForceLogging::logger()){
            step_impl<true>();
        }else{
            step_impl<false>();
        }
    }

    template<bool EnableLogging>
    void step_impl()
    {
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        if(EnableLogging && ForceLogging::logger()){
            ForceLogging::logger()->SetTime(m_state->t);
            ForceLogging::logger()->LogProperty("dt", 1, &m_state->dt);
            double seed_low=m_state->seed &0xFFFFFFFFul;
            double seed_high=m_state->seed>>32;
            ForceLogging::logger()->LogProperty("seed_lo", 1, &seed_low);
            ForceLogging::logger()->LogProperty("seed_high", 1, &seed_high);
            double t_hash_low=m_t_hash&0xFFFFFFFFul;
            double t_hash_high=m_t_hash>>32;
            ForceLogging::logger()->LogProperty("t_hash_lo", 1, &t_hash_low);
            ForceLogging::logger()->LogProperty("t_hash_high", 1, &t_hash_high);
            for(auto &b : m_state->beads){
                double h=b.get_hash_code().hash;
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(), "b_hash", 1, &h);
                double x[3]={b.x[0],b.x[1],b.x[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x",3,x);
                double v[3]={b.v[0],b.v[1],b.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f",3,f);
            }
        }

        double dt=m_state->dt;

        // Clear all cell information
        m_cells.resize(calc_num_cells());
        for(auto &c : m_cells){
            c.clear();
        }

        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &b : m_state->beads){
            if(m_fixed_water){
                if(b.bead_type){
                    dpd_maths_core_half_step::update_pos(dt, m_state->box, b);
                }
            }else{
                dpd_maths_core_half_step::update_pos(dt, m_state->box, b);
            }
            m_cells.at( world_pos_to_cell_index(b.x) ).push_back(&b);
        }

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        vec3i_t pos, dir;
        //std::cerr<<"naive_dpd_engine_half_step : skippping dpd forces.\n";
        
        for(pos[0]=0; pos[0]<m_lengths[0]; pos[0]++){
            for(pos[1]=0; pos[1]<m_lengths[1]; pos[1]++){
                for(pos[2]=0; pos[2]<m_lengths[2]; pos[2]++){
                    auto &home=m_cells[cell_pos_to_index(pos)];
                    for(dir[0]=-1; dir[0]<=+1; dir[0]++){
                        for(dir[1]=-1; dir[1]<=+1; dir[1]++){
                            for(dir[2]=-1; dir[2]<=+1; dir[2]++){
                                auto [other_pos,other_delta] = make_relative_cell_pos(pos, dir);
                                const auto &other=m_cells[cell_pos_to_index(other_pos)];
                                update_cell_forces<EnableLogging>(home, other, other_delta);
                            }
                        }
                    }
                }
            }
        }
        

        // Update all angle bonds
        //std::cerr<<"naive_dpd_engine_half_step_tbb : skippping angle forces.\n";
        
        for(const auto &p : m_state->polymers){
            const auto &pt = m_state->polymer_types.at(p.polymer_type);
            for(const auto &bond : pt.bonds){
                update_hookean_bond(p, pt, bond);
            }
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(p, pt, bond_pair);
            }
        }
        

        // Final momentum update
        for(auto &b : m_state->beads){
            if(m_fixed_water){
                if(b.bead_type){
                    dpd_maths_core_half_step::update_mom(dt, b);
                }
            }else{
                dpd_maths_core_half_step::update_mom(dt, b);
            }
        }

        if(EnableLogging && ForceLogging::logger()){
            for(auto &b : m_state->beads){

                if(ForceLogging::logger()){
                    double v[3]={b.v[0],b.v[1],b.v[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v_next",3,v);
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_next",3,f);
                }
            }
        }

        m_state->t += 1;
    }

    template<bool EnableLogging,class TBeadPtrVec>
    void  __attribute__((noinline))  update_cell_forces(TBeadPtrVec &home, const TBeadPtrVec &other, const vec3r_t &other_delta)
    {
        /*
        for(const Bead *ob : other)
        {
            _mm_prefetch((const char *)&ob->x, _MM_HINT_T0);
        }
        */

        
        for(const Bead *ob : other)
        {
            if(m_fixed_water){
                if(ob->bead_type==0){
                    continue;
                }
            }

            vec3r_t ob_x = vec3r_t(ob->x) + other_delta;
            for(Bead *hb : home){
                vec3r_t dx =  hb->x - ob_x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                assert(hb!=ob);

                double dr=pow_half(dr2);

                const double kappa=0.0, r0=0.5; // Disable hookean bonds
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

                if(EnableLogging){
                    std::pair<const Bead*,const Bead*> pp(hb,ob);
                    if(!m_seen_pairs.insert(pp).second){
                        fprintf(stderr, "  Seen %u %u twice\n", hb->get_hash_code().hash, ob->get_hash_code().hash);
                    }
                }
            }
        }
    }

    void update_hookean_bond(const Polymer &p, const PolymerType &pt, const Bond &b) const
    {
        unsigned head_index=p.bead_ids.at(b.bead_offset_head);
        unsigned tail_index=p.bead_ids.at(b.bead_offset_tail);
        Bead &head_bead=m_state->beads[head_index];
        Bead &tail_bead=m_state->beads[tail_index];

        vec3r_t dx=calc_distance_from_to(head_bead.x, tail_bead.x);
        double dr=dx.l2_norm();

        if(dr>m_max_bond_length){
            throw std::runtime_error("naive_dpd_engine_half_step : an angle bond has snapped.");
        }

        double dr0=b.r0-dr;
        double hookeanForce=b.kappa*dr0;

        vec3r_t f=dx * (hookeanForce/dr);
        head_bead.f -= f;
        tail_bead.f += f;
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
