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

    const int MAX_POLYMER_LENGTH=63;

    struct PolymerInfo
    {
        const PolymerType *polymer_type;
        std::array<Packed*,MAX_POLYMER_LENGTH> bead_storage;
    };

    // One entry per polymer_id, used to find all beads in a polymer
    // This array is very wasteful if:
    // - The bead order is not optimised (i.e. polymers before monomers)
    // - The polymers are quite small compared to MAX_POLYMER_LENGTH
    // It probably doesn't have terrible locality, and we implement as a flat
    // vector to avoid indirecting through a table.
    std::vector<PolymerInfo> m_polymer_info;

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

        #error "Here"

        float dt=m_state->dt;

        // Move the beads, and then assign to cells based on x(t+dt)
        parallel_for_each(m_cells, 32, [&](Cell &c){
            for(int bi=c.packed.size()-1; bi>=0; bi--){
                Packed &b=c.packed[i];
                dpd_maths_core_half_step::update_pos<float>(dt, m_lengthsf, b);

                unsigned index=world_pos_to_cell_index(b->x);
                if(index!=c.index){
                    //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                    // flush to the backing
                    
                    auto &dst_cell=m_cells[index];
                    dst_cell.packed.push_back(b);  // Conflict groups gaurantee this is safe
                    c.packed[bi]=b;
                    c.packed.pop_back();

                    if(!BeadHash{b.id}.is_monomer()){
                        unsigned polymer_id=BeadHash{b.hash}.get_polymer_id();
                        unsigned polymer_offset=BeadHash{b.hash}.get_polymer_offset();

                        assert(m_polymer_info[polymer_id].bead_storage[polymer_offset]==&b;
                        m_polymer_info[polymer_id].bead_storage[polymer_offset]=&dst_cell.packed.back();
                    }
                }
            }
        });

        // At this point each cell will have most beads in c.packed, and might have some in c.packed_incoming

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        parallel_for_each_cell_blocked([&](Cell *c){ process_cell<EnableLogging>(c); } );

        // Update all bonds
        parallel_for_each(m_polymer_info, m_non_monomer_grain, [&](const PolymerInfo *p){
            const auto &pt = 
            for(const auto &bond : pt.bonds){
                update_bond(*p, pt, bond);
            }
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(*p, pt, bond_pair);
            }
        });

        // Final mom
        parallel_for_each(m_cells, 32, [&](Cell &c){
            for(unsigned i=0; i<c.packed.size(); i++){
                dpd_maths_core_half_step::update_mom(m_state->dt, b);
            }
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

    template<bool EnableLogging>
    void process_cell(Cell *c)
    {
        if(c->is_edge){
            process_cell_neighbours<EnableLogging,true>(*c);
        }else{
            process_cell_neighbours<EnableLogging,false>(*c);
        }
        update_intra_forces_packed<EnableLogging>(*c);
    }

    template<bool EnableLogging,bool IsEdge>
    void process_cell_neighbours(Cell &c)
    {
        assert(c.neighbours.size()>=3);

        int i=c.neighbours.size()-1;
        while(i>=0){
            Cell *curr=c.neighbours[i];
            --i;
            update_inter_forces_packed<EnableLogging,IsEdge>(c, *curr);
        }
    }

};

#endif
