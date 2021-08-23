#ifndef naive_dpd_engine_half_step_tbb_v2_hpp
#define naive_dpd_engine_half_step_tbb_v2_hpp

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

class NaiveDPDEngineHalfStepTBBV2
    : public NaiveDPDEngineHalfStep
{
    // This is a 32-byte structure, so we get two per cache line
    struct BeadView
    {
        uint32_t hash;
        uint32_t bead_id;
        vec3f_t x;
        vec3f_t v;

        uint32_t get_hash_code() const
        { return hash; }

        unsigned get_bead_type() const
        { return bead_hash_get_bead_type(hash); }

        unsigned get_polymer_id() const
        { return bead_hash_get_polymer_id(hash); }

        unsigned get_polymer_offset() const
        { return bead_hash_get_polymer_offset(hash); }
    };
    static_assert(sizeof(BeadView)==32);

    struct LockedCell
    {
        unsigned pos_index;
        tbb::concurrent_vector<uint32_t> new_bead_ids;
        std::vector<BeadView> bead_views;
    };

    struct Interaction
    {
        float con;
        float sqrt_diss;
    };

    std::vector<LockedCell> m_locked_cells;
    int m_num_cells = -1;
    float m_inv_sqrt_dt=nan("");
    std::vector<Interaction> m_interactions;

    using range1d_t = tbb::blocked_range<int>;

    void import_cells()
    {
        check_constraints_and_setup();

        m_num_cells=calc_num_cells();
        m_inv_sqrt_dt=recip_pow_half(m_state->dt);
        int cell_reserve_size=4 * m_state->beads.size() / m_num_cells;

        m_interactions.resize(m_numBeadTypes*m_numBeadTypes);
        for(unsigned i=0; i<m_interactions.size(); i++){
            m_interactions[i].con=m_state->interactions[i].conservative;
            m_interactions[i].sqrt_diss=sqrt(m_state->interactions[i].dissipative);
        }

        // Clear all cell information
        m_locked_cells.resize(m_num_cells);
        tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                m_locked_cells[i].bead_views.clear();
                m_locked_cells[i].new_bead_ids.clear();

                m_locked_cells[i].bead_views.reserve(cell_reserve_size);
                m_locked_cells[i].new_bead_ids.reserve(cell_reserve_size);
            }
        });

        // Place the beads in cells, and pre-correct velocity backwards
        double dt=m_state->dt;
        tbb::parallel_for(range1d_t(0, m_state->beads.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &b = m_state->beads[i];
                dpd_maths_core_half_step::update_mom(-dt, b);
                m_locked_cells.at( world_pos_to_cell_index(b.x) ).new_bead_ids.push_back(b.bead_id); // push into concurrent_vector
            }
        });

        // Put them all in the current bead list and update views
       tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c = m_locked_cells[i];
                assert(c.bead_views.size()==0);
                unsigned nbeads=c.new_bead_ids.size();
                c.bead_views.resize(nbeads);
                for(unsigned i=0; i<nbeads; i++){
                    Bead *b=&m_state->beads.at(c.new_bead_ids[i]);
                    c.bead_views[i].bead_id=b->bead_id;
                    c.bead_views[i].hash=b->get_hash_code();
                    c.bead_views[i].x=vec3f_t(b->x);
                    c.bead_views[i].v=vec3f_t(b->v);
                }
                c.new_bead_ids.clear();
            }
        });
    }

    void export_cells()
    {
        float dt=m_state->dt;
        tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c = m_locked_cells[i];
                assert(c.new_bead_ids.empty());
                for(unsigned j=0; j<c.bead_views.size(); j++){
                    auto &b = m_state->beads[c.bead_views[j].bead_id];
                    dpd_maths_core_half_step::update_mom(dt, b);
                }
                c.bead_views.clear();
            }
        });
    }

    void Run(unsigned nSteps) override
    {
        import_cells();

        for(unsigned i=0; i<nSteps; i++){
            step();
        }

        export_cells();
    }

    void step() override
    {
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        float dt=m_state->dt;
        vec3f_t box(m_state->box);

        // Update mom (prev step), then pos (this step), then move if necessary
        // Pre: c.new_beads.empty() for all cells
         tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c =m_locked_cells[i];
                //std::cerr<<"  c.size="<<c.bead_views.size()<<"\n";
                for(int j=c.bead_views.size()-1; j>=0; j--){
                    auto &bv=c.bead_views[j];
                    auto &b=m_state->beads[bv.bead_id];
                    dpd_maths_core_half_step::update_mom(dt, b);
                    //std::cerr<<"  t="<<m_state->t<<", update "<<bv.bead_id<<"\n";
                    dpd_maths_core_half_step::update_pos(dt, box, b);
                    auto pos_index=world_pos_to_cell_index(b.x);
                    if(pos_index!=c.pos_index){
                        m_locked_cells.at(pos_index).new_bead_ids.push_back(bv.bead_id);
                        std::swap(c.bead_views[j], c.bead_views.back());
                        c.bead_views.resize(c.bead_views.size()-1);
                    }
                }
            }
        });

        tbb::parallel_for(range1d_t(0, m_locked_cells.size(), 1024), [&](const range1d_t &r){
            for(int i=r.begin(); i<r.end(); i++){
                auto &c =m_locked_cells[i];
                for(unsigned i=0; i<c.new_bead_ids.size(); i++){
                    Bead &b=m_state->beads[c.new_bead_ids[i]];
                    c.bead_views.resize(c.bead_views.size()+1);
                    c.bead_views.back().bead_id=b.bead_id;
                    c.bead_views.back().hash=b.get_hash_code();
                    c.bead_views.back().x=vec3f_t(b.x);
                    c.bead_views.back().v=vec3f_t(b.v);
                }
                c.new_bead_ids.clear();
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
                                    auto [other_pos,other_delta_d] = make_relative_cell_pos(pos, dir);
                                    auto &other=m_locked_cells[cell_pos_to_index(other_pos)];
                                    update_cell_forces_inter(home.bead_views, other.bead_views, vec3f_t(other_delta_d));
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

        m_state->t += 1;
    }

    bool is_bonded(
        const BeadView &home, const BeadView &other,
        float &kappa,
        float &r0
    ) const
    {
        kappa=0;
        r0=0;
        if(home.get_polymer_id()!=other.get_polymer_id()){
            return false; // Rejects the vast majority
        }
        // This is very slow. In practise use a better method
        Bead &hbb=m_state->beads[home.bead_id];
        for(const Bond &b : m_state->polymer_types.at(hbb.polymer_type).bonds){
            if( ((b.bead_offset_head==home.get_polymer_offset()) && (b.bead_offset_tail==other.get_polymer_offset()))
                ||
                ((b.bead_offset_head==other.get_polymer_offset()) && (b.bead_offset_tail==home.get_polymer_offset()))
            ){
                kappa=b.kappa;
                r0=b.r0;
                return true;
            }
        }
        return false;
    }

    void update_cell_forces_inter(std::vector<BeadView> &home, std::vector<BeadView> &other, const vec3f_t &other_delta) const
    {
        for(BeadView & hb : home){
            for(const BeadView &ob : other)
            {
                vec3f_t dx = hb.x - ob.x - other_delta; 
                float dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1.0f || dr2<0.00001f){
                    continue;
                }

                float dr=pow_half(dr2);

                float kappa,r0;
                is_bonded(hb, ob, kappa, r0);

                vec3f_t f;

                Bead &hbb=m_state->beads[hb.bead_id];
                
                dpd_maths_core_half_step::calc_force<float,vec3f_t>(
                    m_inv_sqrt_dt,
                    [&](unsigned a, unsigned b){ return m_interactions[a*m_numBeadTypes+b].con; },
                    [&](unsigned a, unsigned b){ return m_interactions[a*m_numBeadTypes+b].sqrt_diss; },
                    m_t_hash,
                    dx, dr,
                    kappa, r0,
                    hb,
                    ob,
                    f
                );
                hbb.f += vec3r_t(f);
            }
        }
    }
};

#endif
