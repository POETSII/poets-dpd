#ifndef naive_dpd_engine_half_step_tbb_hpp
#define naive_dpd_engine_half_step_tbb_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_step.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"

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
class NaiveDPDEngineHalfStepTBB
    : public NaiveDPDEngineHalfStep
{
    struct LockedCell
    {
        tbb::concurrent_vector<Bead*> new_beads;
        std::vector<Bead*> beads;
    };

    std::vector<LockedCell> m_locked_cells;

    void step() override
    {
        using range1d_t = tbb::blocked_range<int>;

        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        double dt=m_state->dt;

        int num_cells=calc_num_cells();
        int cell_reserve_size=4 * m_state->beads.size() / num_cells;

        // Clear all cell information
        m_locked_cells.resize(num_cells);
        tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
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

         tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c =m_locked_cells[i];
                c.beads.assign(c.new_beads.begin(), c.new_beads.end());
                c.new_beads.clear();
            }
         });

        // Calculate all the DPD and 2-bead bond forces
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
                                    update_cell_forces(home.beads, other.beads, other_delta);
                                }
                            }
                        }
                    }
                }
            }
        });

        // Update all angle bonds
        // We can do this in parallel, even though update_angle_bond modifies 3 beads, because
        // each polymer is processed in the same task. So there cant be any angle bonds affecting
        // inter-task, as that would mean angle bonds between beads in different polymers.
        tbb::parallel_for(range1d_t(0, m_state->polymers.size(), 64), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                const auto &p = m_state->polymers[i];
                const auto &pt = m_state->polymer_types.at(p.polymer_type);
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

        m_state->t += 1;
    }
};

#endif
