#ifndef naive_dpd_engine_half_merge_tbb_v3_hpp
#define naive_dpd_engine_half_merge_tbb_v3_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_merge.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <immintrin.h>
#include <array>

class NaiveDPDEngineHalfMergeTBBV3
    : public NaiveDPDEngineHalfMerge
{
public:
    virtual double GetMaxBondLength() const
    {
        return 1000;
    }

    void Attach(WorldState *s) override
    {
        NaiveDPDEngineHalfMerge::Attach(s);
        if(s){
            m_lengthsf=m_lengths;
        }
    }

private:
    vec3f_t m_lengthsf;

    struct SuperCell
    {
        std::array<Cell*,8> members;
    };

    virtual void step()
    {
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

        float dt=m_state->dt;

        // Move the beads, and then assign to cells based on x(t+dt)
        parallel_for_each(m_cells, 256, [&](Cell &c){
            for(int bi=c.packed.size()-1; bi>=0; bi--){
                auto &b=c.packed[i];
                dpd_maths_core_half_step::update_pos<float>(dt, m_lengthsf, b);

                //std::cerr<<"In "<<c.pos<<" at "<<m_state->t<<"\n";
                dpd_maths_core_half_step::update_pos(dt, m_lengths, *b);
                unsigned index=world_pos_to_cell_index(b->x);
                if(index!=c.index){
                    //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                    // flush to the backing
                    Bead *backing=c.beads[bi];
                    backing->x=b.x;
                    backing->v=b.v;
                    backing->f=b.f;
                    auto &dst_cell=m_cells[index];
                    {
                        std::unique_lock<std::mutex> lk(dst_cell.mutex);
                        dst_cell.incoming.push_back(b);
                    }
                    c.beads[bi]=c.beads.back();
                    c.packed[bi]=b;
                    c.beads.pop_back();
                    c.packed.pop_back();
                }
            }
        });

        // At this point each cell will have most beads in c.beads, and might have some in c.incoming

        // Turns out more efficient to do explicitly than all at once
        parallel_for_each(m_cells, 2048, [&](Cell &c){
            transfer_incoming(c);
        });

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        parallel_for_each_cell_blocked([&](Cell *c){ process_cell<EnableLogging>(c); } );

        // Update all bonds
        if(0){
            parallel_for_each(m_non_monomers, m_non_monomer_grain, [&](const Polymer *p){
                const auto &pt = m_state->polymer_types.at(p->polymer_type);
                for(const auto &bond : pt.bonds){
                    update_bond(*p, pt, bond);
                }
                for(const auto &bond_pair : pt.bond_pairs){
                    update_angle_bond(*p, pt, bond_pair);
                }
            });
        }else{
            for(const auto &p : m_state->polymers){
                const auto &pt = m_state->polymer_types.at(p.polymer_type);
                for(const auto &bond : pt.bonds){
                    update_bond(p, pt, bond);
                }
                for(const auto &bond_pair : pt.bond_pairs){
                    update_angle_bond(p, pt, bond_pair);
                }
            }
        }

        // Final mom
        parallel_for_each(m_state->beads, 1024, [&](Bead &b){
            dpd_maths_core_half_step::update_mom(m_state->dt, b);
        });

        if(EnableLogging && ForceLogging::logger()){
            for(auto &b : m_state->beads){

                if(ForceLogging::logger()){
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_next",3,f);
                }
            }
        }

        m_state->t += 1;
    }

    // Idempotent function to moving incoming to beads. Should be fast in case where incoming is empty
    void transfer_incoming(Cell &c)
    {
        if(!c.incoming.empty()){
            c.beads.insert(c.beads.end(), c.incoming.begin(), c.incoming.end());
            c.incoming.clear();
        }
    };

    template<bool EnableLogging>
    void process_cell(Cell *c)
    {
        if(c->is_edge){
            process_cell_neighbours<EnableLogging,true>(*c);
        }else{
            process_cell_neighbours<EnableLogging,false>(*c);
        }
        update_intra_forces<EnableLogging>(*c);
    }

    template<bool EnableLogging,bool IsEdge>
    void process_cell_neighbours(Cell &c)
    {
        assert(c.neighbours.size()>=3);

        int i=c.neighbours.size()-1;
        Cell *next_next_next=c.neighbours[i-2];
        Cell *next_next=c.neighbours[i-1];
        Cell *next=c.neighbours[i];
        while(i>=0){
            Cell *curr=next;
            next=next_next;
            next_next=next_next_next;
            //Cell *curr=c.neighbours[i];
            --i;
            if(i>=2){
                next_next_next=c.neighbours[i-2];
                _mm_prefetch(&next_next_next->beads, _MM_HINT_T0);
            }
            if(i>=1){
                _mm_prefetch(&next_next->beads[0], _MM_HINT_T0);
            }
            if(next){
                for(auto p : next->beads){
                    _mm_prefetch(&p->x, _MM_HINT_T0);
                }
            }
            update_inter_forces<EnableLogging,IsEdge>(c, *curr);
        }
    }

};

#endif
