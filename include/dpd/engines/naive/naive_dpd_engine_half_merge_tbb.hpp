#ifndef naive_dpd_engine_half_merge_tbb_hpp
#define naive_dpd_engine_half_merge_tbb_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_merge.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <immintrin.h>
#include <array>

class NaiveDPDEngineHalfMergeTBB
    : public NaiveDPDEngineHalfMerge
{
public:
    void Attach(WorldState *s) override
    {
        NaiveDPDEngineHalfMerge::Attach(s);
        if(s){
            make_conflict_groups();
            collect_non_monomers();
        }
    }

private:
    /*
    struct MiniCell
    {
        uint32_t resident_beads;
        __m128 x_packed[3];
        __m128 v_packed[3];
        __m128 f_packed[3];
        __m128 id_packed;

        std::vector<Bead*> beads;

        void interact(float x[3], float v[3], uint32_t id, float f[3])
        {
            __m128 x_x=_mm_set1_ps(x[0]);
            __m128 x_y=_mm_set1_ps(x[1]);
            __m128 x_z=_mm_set1_ps(x[2]);

            __m128 dx_x=_mm_sub_ps(x_packed[0], x_x);
            __m128 dx_y=_mm_sub_ps(x_packed[1], x_y);
            __m128 dx_z=_mm_sub_ps(x_packed[2], x_z);

            __m128 dr = _mm_add_ps(dx_x, _mm_add_ps(dx_y, dx_z));
            __m128 valid = _mm_and_ps( _mm_cmple_ps( dr, _mm_set1_ps(UPPER_R) ) , _mm_cmpge_ps( dr, _mm_set1_ps(LOWER_R)) );
            
            int mask=_mm_movemask_ps(valid);
            if(mask==0){
                return;
            }
        }
    };
    */

    struct SuperCell
    {
        std::array<Cell*,8> members;
    };

    std::vector<std::vector<SuperCell>> m_conflict_groups;

    std::vector<Polymer*> m_non_monomers;
    size_t m_non_monomer_grain;

    void collect_non_monomers()
    {
        m_non_monomers.clear();
        size_t total_bonds=0, total_bond_pairs=0;
        for(Polymer &p : m_state->polymers){
            if(p.bead_ids.size()>1){
                m_non_monomers.push_back(&p);
                PolymerType &pt=m_state->polymer_types.at(p.polymer_type);
                total_bonds += pt.bonds.size();
                total_bond_pairs += pt.bond_pairs.size();
            }
        }
        if(m_non_monomers.empty()){
            m_non_monomer_grain=0;
        }else{
            double avg_ops_per_polymer=(total_bonds * 30 + total_bond_pairs * 50)/m_non_monomers.size();
            m_non_monomer_grain=(unsigned)std::max(1.0, 1000000 / avg_ops_per_polymer);
        }
    }

    /*void make_conflict_groups()
    {
        assert(!m_cells.empty());
        m_conflict_groups.clear();

        for(int ox=0; ox<4; ox++){
            for(int oy=0; oy<4; oy++){
                for(int oz=0; oz<4; oz++){
                    m_conflict_groups.push_back({});
                    auto &group=m_conflict_groups.back();

                    for(int ix=ox; ix<m_dims[0] ; ix+=4){
                        for(int iy=oy; iy<m_dims[1] ; iy+=4){
                            for(int iz=oz; iz<m_dims[2] ; iz+=4){
                                group.push_back( &m_cells.at(cell_pos_to_index({ix,iy,iz})));
                            }
                        }
                    }
                }
            }
        }
    }*/

    void make_conflict_groups()
    {
        std::unordered_set<unsigned> seen;

        m_conflict_groups.clear();
        for(unsigned gx=0; gx<4; gx+=2){
            for(unsigned gy=0; gy<4; gy+=2){
                for(unsigned gz=0; gz<4; gz+=2){
                    // This gives the origin of a 2x2x2 cube within a 4x4x4 block
                    std::vector<SuperCell> group;
                    // Loop over all super cells within group
                    for(unsigned ix=gx; ix<m_dims[0]; ix+=4){
                        for(unsigned iy=gy; iy<m_dims[1]; iy+=4){
                            for(unsigned iz=gz; iz<m_dims[2]; iz+=4){
                                SuperCell sc;
                                unsigned off=0;
                                for(int lx=0; lx<2; lx++){
                                    for(int ly=0; ly<2; ly++){
                                        for(int lz=0; lz<2; lz++){
                                            unsigned index=cell_pos_to_index({ix+lx,iy+ly,iz+lz});
                                            if(!seen.insert(index).second){
                                                throw std::runtime_error("Duplicate");
                                            }
                                            sc.members[off++]=&m_cells[index];
                                        }
                                    }
                                }
                                group.push_back(sc);
                            }
                        }
                    }
                    m_conflict_groups.push_back(group);
                }
            }
        }
    }

    virtual void step()
    {
        if(ForceLogging::logger()){
            step_impl<true>();
        }else{
            step_impl<false>();
        }
    }

    template<class T,class F>
    void parallel_for_each(std::vector<T> &x, unsigned grain, F &&f)
    {
        using range_t=tbb::blocked_range<size_t>;
        tbb::parallel_for(range_t(0,x.size(),grain), [&](const range_t &r){
            for(size_t i=r.begin(); i<r.end(); i++){
                f(x[i]);
            }
        }, tbb::simple_partitioner{});
    }

    /*template<class F>
    void parallel_for_each_cell_blocked(unsigned grain, F &&f)
    {
        for(unsigned i=0; i<m_conflict_groups.size(); i++){
            parallel_for_each(m_conflict_groups[i], grain, f);
        }
    }*/

    template<class F>
    void parallel_for_each_cell_blocked(F &&f)
    {
        for(unsigned i=0; i<m_conflict_groups.size(); i++){
            parallel_for_each(m_conflict_groups[i], 1, [&](const SuperCell &c) {
                for(unsigned i=0; i<8; i++){
                    f(c.members[i]);
                }
            });
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
                Bead bt=b;
                dpd_maths_core_half_step::update_mom(m_state->dt/2, bt);
                double v[3]={bt.v[0],bt.v[1],bt.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f",3,f);
            }
        }

        double dt=m_state->dt;

        // Move the beads, and then assign to cells based on x(t+dt)
        parallel_for_each(m_cells, 256, [&](Cell &c){
            for(int bi=c.beads.size()-1; bi>=0; bi--){
                Bead *b=c.beads[bi];
                //std::cerr<<"In "<<c.pos<<" at "<<m_state->t<<"\n";
                dpd_maths_core_half_step::update_pos(dt, m_lengths, *b);
                unsigned index=world_pos_to_cell_index(b->x);
                if(index!=c.index){
                    //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                    auto &dst_cell=m_cells[index];
                    {
                        std::unique_lock<std::mutex> lk(dst_cell.mutex);
                        dst_cell.incoming.push_back(b);
                    }
                    c.beads[bi]=c.beads.back();
                    c.beads.pop_back();
                }
            }
        });

        // At this point each cell will have most beads in c.beads, and might have some in c.incoming

        // Turns out more efficient to do explicitly than all at once
        parallel_for_each(m_cells, 1024, [&](Cell &c){
            transfer_incoming(c);
        });

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        parallel_for_each_cell_blocked([&](Cell *c){ process_cell<EnableLogging>(c); } );

        // Update all bonds
        parallel_for_each(m_non_monomers, m_non_monomer_grain, [&](const Polymer *p){
            const auto &pt = m_state->polymer_types.at(p->polymer_type);
            for(const auto &bond : pt.bonds){
                update_bond(*p, pt, bond);
            }
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(*p, pt, bond_pair);
            }
        });

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
