#ifndef tolerant_dpd_engine_half_step_tbb_hpp
#define tolerant_dpd_engine_half_step_tbb_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_step.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"
#include "dpd/core/logging.hpp"

#include <cassert>
#include <cmath>
#include <array>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range3d.h"
#include "tbb/concurrent_vector.h"

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
class TolerantDPDEngineHalfStepTBB
    : public NaiveDPDEngineHalfStep
{
    struct LockedCell
    {
        unsigned index;
        tbb::concurrent_vector<Bead*> new_beads;
        std::vector<Bead*> beads;
    };

    std::vector<LockedCell> m_locked_cells;

    std::vector<Polymer*> m_polymers_with_bonds;

public:
    virtual double GetMaxBondLength() const
    {
        return 1000;
    }

    void Attach(WorldState *s) override
    {
        m_locked_cells.clear();
        NaiveDPDEngineHalfStep::Attach(s);

        m_polymers_with_bonds.clear();
        if(s){
            for(auto &p : s->polymers){
                auto &pt=s->polymer_types.at(p.polymer_type);
                if(!pt.bonds.empty()){
                    m_polymers_with_bonds.push_back(&p);
                }
            }
        }
    }
protected:

    template<bool EnableLogging>
    void update_hookean_bond(const Polymer &p, const PolymerType &pt, const Bond &b) const
    {
        unsigned head_index=p.bead_ids[b.bead_offset_head];
        unsigned tail_index=p.bead_ids[b.bead_offset_tail];
        Bead &head=m_state->beads[head_index];
        Bead &tail=m_state->beads[tail_index];

        vec3r_t dx=calc_distance_from_to(head.x, tail.x);
        double r=dx.l2_norm();

        vec3r_t head_f, tail_f;

        dpd_maths_core_half_step::calc_hookean_force<EnableLogging>(
            b.kappa, b.r0, 
            dx, r, 
            head_f, tail_f
        );
        
        head.f += head_f;
        tail.f += tail_f;
        
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
            vec3r_t ob_x = vec3r_t(ob->x) + other_delta;
            for(Bead *hb : home){
                vec3r_t dx =  hb->x - ob_x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                assert(hb!=ob);

                double dr=pow_half(dr2);

                double kappa=0.0, r0=0.5;

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

    void step() override
    {
        if(ForceLogging::logger()){
            ForceLogging::logger()->SetTime(m_state->t);
        }

        using range1d_t = tbb::blocked_range<int>;

        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        double dt=m_state->dt;

        int num_cells=calc_num_cells();

        if(m_locked_cells.empty()){
            int cell_reserve_size=4 * m_state->beads.size() / num_cells;

            // Clear all cell information
            m_locked_cells.resize(num_cells);
            tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 64), [&](const range1d_t &r){
                for(int i=r.begin(); i<r.end(); i++){
                    m_locked_cells[i].index=i;
                    m_locked_cells[i].beads.clear();
                    m_locked_cells[i].beads.reserve(cell_reserve_size);
                }
            });

            // Move the beads, and then assign to cells based on x(t+dt)
            tbb::parallel_for(range1d_t(0, m_state->beads.size(), 1024), [&](const range1d_t &r){
                for(int i=r.begin(); i<r.end(); i++){
                    auto &b = m_state->beads[i];
                    dpd_maths_core_half_step::update_pos(dt, m_state->box, b);
                    m_locked_cells.at( world_pos_to_cell_index(b.x) ).new_beads.push_back(&b); // push into concurrent_vector
                }
            });
        }else{
            // Move the beads, and then assign to cells based on x(t+dt)
            tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 64), [&](const range1d_t &r){
                for(int i=r.begin(); i<r.end(); i++){
                    auto &c = m_locked_cells[i];
                    for(int j=c.beads.size()-1; j>=0; j--){
                        auto *b = c.beads[j];
                        dpd_maths_core_half_step::update_pos(dt, m_state->box, *b);
                        unsigned index=world_pos_to_cell_index(b->x);
                        if(index!=(unsigned)i){
                            // Add to dst
                            assert(m_locked_cells.at(index).index==index);
                            m_locked_cells.at( index ).new_beads.push_back(b); // push into concurrent_vector
                            // Remove from here
                            c.beads[j]=c.beads.back();
                            c.beads.resize(c.beads.size()-1);
                        }else{
                            assert(c.index==index);
                        }
                    }
                }
            });
        }

         tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 64), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c =m_locked_cells[i];
                if(!c.new_beads.empty()){
                    c.beads.insert(c.beads.end(), c.new_beads.begin(), c.new_beads.end());
                    c.new_beads.clear();
                }
            }
         });

        // Calculate all the DPD forces, but _not_ the 2-bead bond forces
        // Each cell's force is calculated independently
        using range3d_t = tbb::blocked_range3d<int,int,int>;
        const int GRAIN_CELL=2;
        range3d_t cell_range(0, m_lengths[0], GRAIN_CELL, 0, m_lengths[1], GRAIN_CELL, 0, m_lengths[2], GRAIN_CELL);
        tbb::parallel_for(cell_range, [&](const range3d_t &r){
            vec3i_t pos, dir;
            for(pos[0]=r.pages().begin(); pos[0]<r.pages().end(); pos[0]++){
                for(pos[1]=r.rows().begin(); pos[1]<r.rows().end(); pos[1]++){
                    for(pos[2]=r.cols().begin(); pos[2]<r.cols().end(); pos[2]++){
                        auto &home=m_locked_cells[cell_pos_to_index(pos)];
                        for(dir[0]=-1; dir[0]<=+1; dir[0]++){
                            for(dir[1]=-1; dir[1]<=+1; dir[1]++){
                                for(dir[2]=-1; dir[2]<=+1; dir[2]++){
                                    auto [other_pos,other_delta] = make_relative_cell_pos(pos, dir);
                                    const auto &other=m_locked_cells[cell_pos_to_index(other_pos)];
                                    update_cell_forces<false>(home.beads, other.beads, other_delta);
                                }
                            }
                        }
                    }
                }
            }
        });

        // Update all hookean and angle bonds
        // We can do this in parallel, even though update_angle_bond modifies 3 beads, because
        // each polymer is processed in the same task. So there cant be any angle bonds affecting
        // inter-task, as that would mean angle bonds between beads in different polymers.
        tbb::parallel_for(range1d_t(0, m_polymers_with_bonds.size(), 64), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                const auto &p = *m_polymers_with_bonds[i];
                const auto &pt = m_state->polymer_types.at(p.polymer_type);
                // TODO : This could be more efficient, combining hookean and angle... But so could verything else
                for(const auto &bond : pt.bonds){
                    update_hookean_bond<false>(p, pt, bond);
                }
                for(const auto &bond_pair : pt.bond_pairs){
                    update_angle_bond(p, pt, bond_pair);
                }
            }
        });

        // Final momentum update
        tbb::parallel_for(range1d_t(0, m_state->beads.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                dpd_maths_core_half_step::update_mom(dt, m_state->beads[i]);
            }
        });

        if(ForceLogging::logger()){
            for(auto &b : m_state->beads){
                double x[3]={b.x[0],b.x[1],b.x[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x_next",3,x);
            }
        }

        m_state->t += 1;
    }
};

#endif
