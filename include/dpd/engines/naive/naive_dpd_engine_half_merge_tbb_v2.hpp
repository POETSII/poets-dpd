#ifndef naive_dpd_engine_half_merge_tbb_hpp
#define naive_dpd_engine_half_merge_tbb_hpp

#error "Not complete or tested"

#include "dpd/engines/naive/naive_dpd_engine_half_merge.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <immintrin.h>
#include <array>

#define SIMDPP_ARCH_X86_AVX2
#include "simdpp/simd.h"

class NaiveDPDEngineHalfMergeTBBV2
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
    using vec4f = simdpp::float32<4>;

    struct vec4f_plus_u32
    {
        union{
            float x[4]; // Stored in 1,2, and 3
            struct {
                uint32_t _pad_[3];
                uint32_t id;
            };
        };

        uint32_t get_u32() const
        { return id; }

        float get_f32(unsigned i) const
        {
            return x[i-1];
        }

        vec4f get_vec4f() const
        {
            vec4f_t res=load(x);
            return simdpp::move4_l<1>(res); // Shift each element towards 0; MSW is zero
        }
    };

    struct LocalBead
    {
        vec3_plus_u32 v_plus_id;
        vec4f f;
    };

    struct Cell
    {

        static constexpr int DIRECT_STORE=4;

        uint16_t resident;
        uint16_t allocated;

        vec4f *x;
        LocalBead *bead_extra;
    };

    void interact_direct(Cell &home, Cell &other)
    {
        vec4f dx_i[MAX_BEADS_PER_CELL];
        float dr2_i[MAX_BEADS_PER_CELL];

        for(int io=0; io < other.resident; io++){
            vec4f other_x=other.x[io];
            vec4f other_f(0,0,0,0);

            // TODO : move after for low probability interactions?
            vec3_plus_u32 other_v_plus_id=other.bead_extra[io].v_plus_id

            bool any_hit=false;
            for(int i=0; i<home.resident; i++){
                vec4f dx = x[i] - other_x;
                vec4f dx2=dx[i]*dx[i];
                float dr2=simdpp::reduce_add(dx2);
                dx_i[i] = dx;
                dr2_i[i] = dr2;
                any_hit |= dr2 < 1;
            }

            if(!any_hit){
                return;
            }

            BeadHash other_bead_hash{other_v_plus_id.get_u32()};

            const float *other_strengths_row = m_interactions + 2*MAX_BEAD_TYPES*other_bead_hash.get_bead_type();

            for(int i=0; i<home.resident; i++){
                float dr2 = dr2_i[i];
                if(dr2[i] >= 1){
                    continue;
                }
                if(dr2[i]<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                vec4f dx=dx_i[i];

                float dr=sqrtf(dr2);
                float inv_dr=1.0f/dr;

                auto &home_bead_extra=home.bead_extra[i];
                vec3_plus_u32 home_v_plus_id=home_bead_extra.v_plus_id;
                BeadHash home_bead_hash{home_v_plus_id.get_u32()};

                const float *strengths_cell = other_strengths_row + 2*home_bead_hash.get_bead_type();

                vec4f dv = home_v_plus_id - other_v_plus_id;

                float wr = 1.0f - dr;

                float conForce=strengths_cell[0] * wr;

                dx *= inv_dr;

                float rdotv = libsimdpp::reduce_add(dx * dv);
                float sqrt_gammap = strengths_cell[1] * wr;

                float dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
                float u = dpd_maths_core::default_hash(t_hash, home_hash, other_hash);
                float randScale = sqrt_gammap * m_scale_inv_sqrt_dt ;
                float randForce = randScale * u;

                float scaled_force = conForce + dissForce + randForce;

                vec4f f = dx * scaled_force;

                home_bead_extra.f += f;

                other_f -= f;
            }
        }
    }

    

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
