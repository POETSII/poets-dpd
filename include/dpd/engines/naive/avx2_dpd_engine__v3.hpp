


#define VCL_NAMESPACE vcl
#include "dpd/external/vcl/vectorclass.h"

#include <vector>
#include <type_traits>
#include <array>


#include "dpd/maths/dpd_maths_core_half_step.hpp"
#include "dpd/core/dpd_state.hpp"

class AVX2DPDEngineV3
{
    struct w4{
        static const unsigned W = 4;
        using float32 = vcl::Vec4f;
        using uint32= vcl::Vec4ui;
        using uint64 = vcl::Vec4uq;
    };

    struct w8{
        static const unsigned W = 8;
        using float32 = vcl::Vec8f;
        using uint32= vcl::Vec8ui;
        using uint64 = vcl::Vec8uq;
    };

    template<unsigned N>
    struct wn{
        static const unsigned W = N;
        static_assert(N<=8);
        using float32 = std::conditional<N<=4,w4::float32,w8::float32>;
        using uint32 = std::conditional<N<=4,w4::uint32,w8::uint32>;
        using uint64 = std::conditional<N<=4,w4::uint32,w8::uint64>;
    };

    static const unsigned MAX_BEAD_TYPES=8;

    struct Packed
    {
        float x[3];
        float v[3];
        float f[3];
        uint32_t hash;
    };

    struct Cell
    {
        unsigned n;
        union{
            struct{
                w4::float32 x[3];
                w4::float32 v[3];
                w4::float32 f[3];
                w4::uint32 hash;
            } len4;
            struct len8{
                w8::float32 x[3];
                w8::float32 v[3];
                w8::float32 f[3];
                w8::uint32 hash;
            } len8;
        };

        template<class TH>
        auto ref_x(unsigned i) -> typename TH::float32&
        { if const(TH::W==4){ return len4.x[i]; } else { return len8.x[i]; }; }

        template<class TH>
        auto ref_x() -> typename TH::float32 *
        { if const(TH::W==4){ return len4.x; } else { return len8.x; }; }

        template<class TH>
        auto ref_v(unsigned i) -> typename TH::float32&
        { if const(W==4){ return len4.v[i]; } else { return len8.v[i]; }; }

        template<class TH>
        auto ref_v() -> typename TH::float32 *
        { if const(W==4){ return len4.v[i]; } else { return len8.v[i]; }; }

        template<class TH>
        auto ref_f(unsigned i) -> typename TH::float32&
        { if const(W==4){ return len4.f[i]; } else { return len8.f[i]; }; }

        template<class TH>
        auto ref_hash() -> typename TH::uint32&
        { if const(W==4){ return len4.hash[i]; } else { return len8.hash[i]; }; }


        struct NeighbourInfo{
            uint32_t cell_index;        // Index into array of cells
            float offset[3];            // Offsets to apply for positions in neighbour
        };

        std::array<NeighbourInfo,13> forward;
        std::vector<Packed> overspill;
    };

    
    WorldState *m_state;
    float m_scale_inv_sqrt_dt;
    uint64_t m_t_hash;
    std::array<float,MAX_BEAD_TYPES*MAX_BEAD_TYPES> m_con_strength;
    std::array<float,MAX_BEAD_TYPES*MAX_BEAD_TYPES> m_sqrt_diss_strength;

    std::vector<Cell> m_cells;
    std::vector<std::vector<std::vector<uint32_t>>> m_conflict_groups; // Offset into cells.

    template<class F>
    void parallel_for_each_cluster(F &&f)
    {
        for(const auto &cg : m_conflict_groups){
            tbb::parallel_for_each(cg, [&](const std::vector<uint32_t> &cluster){
                for(auto index : cluster){
                    f(m_cells[index]);
                }
            });
        }
    }

    void step()
    {
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        parallel_for_each_cluster([&](Cell &){
            
        });
    }

    void interact(Cell &home)
    {
        if(home.n<=4){
            interact_h<w4>(home);
        }else if(home.n<=4){
            interact_h<w8>(home);
        }else{
            throw std::runtime_error("TODO");
        }
    }

    template<class TH>
    void interact_h(Cell &home)
    {
        for(unsigned i=0; i<13; i++){
            const NeighbourInfo &ni=home.forward[i];
            Cell &other=*m_cells[hi.index];
            if(other.n<=4){
                interact_h_o<TH,w4>(home, other);
            }else if(other.n<=8){
                interact_h_o<TH,w8>(home, other);
            }else{
                throw std::runtime_error("TODO");
            }
        }
    }

    template<class THOut, class THIn>
    typename THOut::float32 broadcast_ith(const typename THIn::float32 &src, unsigned i)
    {
        assert(i < THIn::W);
        return THOut::float32( src[i] );
    }

    template<class TH, class TO>
    void interact_h_o(
        Cell &home,
        Cell &other,
        const float other_offset[3]
    ){
        bool hit_o=0;
        for(unsigned i=0; i<home.n; i++){
            TH::float32 dx[3]={
                broadcast_ith(home.ref_x<TH>(0),i) - other.ref_x<TO>(0) + other_offset[0],
                broadcast_ith(home.ref_x<TH>(1),i) - other.ref_x<TO>(1) + other_offset[1],
                broadcast_ith(home.ref_x<TH>(2),i) - other.ref_x<TO>(2) + other_offset[2]
            };
            auto dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            float min_dr2=vcl::horizontal_min(dr2);
            if(min_dr2 >= 1.0f){
                return;
            }

            hit_o = true;

            uint32_t hhash_i = hhash[i];
            unsigned hbead_type = BeadHash{hhash_i}.get_bead_type();
            float hf_i[3];
            void calc_dpd_forces<TO>(
                m_scale_inv_sqrt_dt, m_t_hash,

                hhash_i,
                hv[0][i],
                hv[1][i],
                hv[2][i],

                &m_con_strength[ hbead_type*MAX_BEAD_TYPES ],
                &m_sqrt_diss_strength[ hbead_type*MAX_BEAD_TYPES ],

                other.ref_hash<TO>(),
                other.get_v<TO>(0), 
                other.get_v<TO>(1), 
                other.get_v<TO>(2), 

                vcl::sqrt(dr2),
                dx[3],

                hf_i,              // in/out : updated in place
                other.ref_f<TO>()
            );
            for(int d=0; d<3; d++){
                hf[d][i] += hf_i[d]; // TODO: convert to vectors, or will compiler do it?
            }
        }
    }

    

    template<class TO>
    static void calc_dpd_forces(
        float scale_inv_sqrt_dt,
        const uint64_t t_hash,

        const uint32_t hhash,
        const float hv_0,
        const float hv_1,
        const float hv_2,

        // These are all indexed by the bead type of home bead
        const float *conStrengthRow,
        const float *sqrtDissStrengthRow,

        const typename TO::uint32_t ohash,
        const typename TO::float32 ov_0,
        const typename TO::float32 ov_1,
        const typename TO::float32 ov_2,

        const typename TO::float32 dr,
        typename TO::float32 dx[3],

        float hf[3],              // in/out : updated in place
        typename TO::float32 of[3] // out: always overwritten
    ){
        assert( vcl::horizontal_and(home_hash != other_hash) );
        assert( vcl::horizontal_and( 0 < dr && dr < 1) );

        auto inv_dr = 1.0f/dr;
            
        typename TO::float32 dv[3]={ hv_0-ov_0, hv_1-ov_1, hv_2-ov_2 };
            
        auto wr = (1.0f - dr);

        auto beadTypes=BeadHash::get_bead_type(ohash);
        auto conStrength = vcl::lookup<MAX_BEAD_TYPES>(beadTypes,conStrengthRow);
        auto sqrtDissStrength = vcl::lookup<MAX_BEAD_TYPES>(beadTypes,sqrtDissStrengthRow);
            
        auto conForce = conStrength*wr;

        for(int d=0; d<3; d++){
            dx[d] = dx[d] * inv_dr;
        }
            
        auto rdotv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
        auto sqrt_gammap = sqrtDissStrength*wr;

        auto dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
        auto u = dpd_maths_core::default_hash(t_hash, home_hash, other_hash);
        auto randForce = sqrt_gammap * scale_inv_sqrt_dt * u;

        auto scaled_force = conForce + dissForce + randForce;

        for(int d=0; d<3; d++){
            auto f = dx[d] * scaled_force;
            hf[d] += vcl::horizontal_add(f);
            of[d] -= f;
        }
    }

};