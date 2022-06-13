#ifndef dpd_maths_core_simd_hpp
#define dpd_maths_core_simd_hpp

#include <immintrin.h>

#include "dpd/core/hash.hpp"
#include "dpd/maths/dpd_maths_core.hpp"

#ifndef PDPD_TINSEL
#include <iostream>
#endif

namespace dpd_maths_core_simd
{

    using dpd_maths_core::default_hash;

    /*  This uses a 64-bit XorShift, with additive combination
        of previous and next values. In scalar terms it is:

        uint64_t xn = XorShift64(xp);
        uint32_t a = (xn>>32) + (xp&0xFFFFFFFF);
        uint32_t b = (xp>>32) + (xn&0xFFFFFFFF);

        This is "good enough". It passes SmallCrush both horizontally
        and vertically, and has a 64-bit period. Full 256-bit
        would be better, but this works well enough and is
        short.

        A mildly suspicious value (p=0.9994) in horizontal Crush todo with weight
        distribution, but nothing close to catastrophic.
        There is 1 catastrophic failure in vertical Crush for linear complexity
        in bit 29, but that is expected given the underlying linear
        generator and the additive mixing. Apart from that it is fine.

        Instructions in avx2 are:
        f(long long __vector(4)&):                             # @f(long long __vector(4)&)
            vmovdqa ymm0, ymmword ptr [rdi]
            vpsllq  ymm1, ymm0, 13
            vpxor   ymm1, ymm1, ymm0
            vpsrlq  ymm2, ymm1, 7
            vpxor   ymm1, ymm2, ymm1
            vpsllq  ymm2, ymm1, 17
            vpxor   ymm1, ymm2, ymm1
            vmovdqa ymmword ptr [rdi], ymm1
            vpshufd ymm1, ymm1, 177                 # ymm1 = ymm1[1,0,3,2,5,4,7,6]
            vpaddd  ymm0, ymm1, ymm0
            ret

        Load and store can often be optimised out, so it is about 8 instructions,
        or 1 instruction per 32-bit number.
    */
    static __m256i XorShift64Add(__m256i &x)
    {
        __m256i a=x;

        x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
        x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
        x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));

        __m256i b=_mm256_shuffle_epi32(x, (2<<6)|(3<<4)|(0<<2)|(1<<0));

        /*x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
        x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
        x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));    

        __m256i c=x;*/

        return _mm256_add_epi32(a,b);
    }

    // TODO : Proper mixing version
    static __m256 XorShift32(__m256i &state)
    {
        __m256i init=state;

        state = state ^ _mm256_slli_epi32( state , 13 );
        state = state ^ _mm256_srli_epi32( state , 7 );
        state = state ^ _mm256_slli_epi32( state , 5 );

        __m256i x = _mm256_add_epi32(init , state);
        const __m256 scale = _mm256_set1_ps(0.00000000023283064365386962890625f);
        __m256 xf = _mm256_cvtepi32_ps(x);
        return xf*scale;
    }

    // TODO : Be a bit critical...
    /*static __m256 XorShift32Weyl(__m256i state[2])
    {
        __m256i init=state;

        state[0] = state[0] ^ _mm256_slli_epi32( state[0] , 13 );
        state[0] = state[0] ^ _mm256_srli_epi32( state[0] , 7 );
        state[0] = state[0] ^ _mm256_slli_epi32( state[0] , 5 );

        // This is actually a quasi-random weyl generator, so is overly even...
        static const uint32_t step=uint32_t(2654435769); // = (1.0/1.6180339887498948482) * 0xFFFFFFFFul;
        state[1] = state[1] + _mm256_set1_epi32(step)

        __m256i x = (init + state[0]) ^ state[1];
        const __m256 scale = _mm256_set1_ps(0.00000000023283064365386962890625f);
        __m256 xf = _mm256_cvtepi32_ps(x);
        return xf*scale;
    }*/

    /*static __m256 XorShift64(__m256i state[2])
    {
        __m256i init=state;

        state[0] = state[0] ^ _mm256_slli_epi64( state , 13 );
        state[0] = state[0] ^ _mm256_srli_epi64( state , 7 );
        state[0] = state[0] ^ _mm256_slli_epi64( state , 17 );
        state[1] = state[1] ^ _mm256_slli_epi64( state , 13 );
        state[1] = state[1] ^ _mm256_srli_epi64( state , 7 );
        state[1] = state[1] ^ _mm256_slli_epi64( state , 17 );

        __m256i x = _mm256_add_epi32(state[0], state[1]);
        const __m256 scale = _mm256_set1_ps(0.00000000023283064365386962890625f);
        __m256 xf = _mm256_cvtepi32_ps(x);
        return xf*scale;
    }*/

    // https://stackoverflow.com/a/23190168
    // TODO : not nesc. the fastest with newer architectures
    static inline float _mm256_reduce_add_ps(__m256 x) {
        /* ( x3+x7, x2+x6, x1+x5, x0+x4 ) */
        const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(x, 1), _mm256_castps256_ps128(x));
        /* ( -, -, x1+x3+x5+x7, x0+x2+x4+x6 ) */
        const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
        /* ( -, -, -, x0+x1+x2+x3+x4+x5+x6+x7 ) */
        const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
        /* Conversion to float is a no-op on x86-64 */
        return _mm_cvtss_f32(x32);
    }

    /* Returns true if there is any interaction (non-zero force), false otherwise.
        The f values are only valid iff it returns true. f is not set to zero
        if there are no interactions.

        This function treats all home lanes as active. It is up to the caller to
        make sure they are valid (e.g. by moving inactive lanes to an x value which
        can't interact) or by applying masking externally.
    */
    template<unsigned MAX_BEAD_TYPES>
    static bool interact_vec8_to_scalar(
        const float scale_inv_sqrt_dt,
        const float conservative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        const float sqrt_dissipative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        
        __m256i &rng_state,
        
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],

        unsigned other_bead_type,
        float other_x[3],
        float other_v[3],

        __m256 f[3] // Force on home beads
    ){
        static_assert(MAX_BEAD_TYPES <= 8);

        __m256 dx[3];
        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            dx[d] = home_x[d] - _mm256_set1_ps(other_x[d]);
        }
        __m256 dr2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        __m256i active = (dr2 > 0.00001f && 1.0f > dr2);
        if(0==_mm256_movemask_epi8(active)){
           return false;
        }

        __m256 dr = _mm256_sqrt_ps (dr2);
        __m256 inv_dr = 1.0f / dr;

        __m256 con_strength_row = _mm256_loadu_ps( conservative_matrix + 8 * other_bead_type);
        __m256 sqrt_diss_strength_row = _mm256_loadu_ps( sqrt_dissipative_matrix + 8 * other_bead_type);
        __m256 con_strength = _mm256_permutevar8x32_ps( con_strength_row, home_bead_types );
        __m256 sqrt_diss_strength = _mm256_permutevar8x32_ps( sqrt_diss_strength_row, home_bead_types );

        auto wr = 1.0f - dr;
        auto conForce = con_strength * wr;

        __m256 dv[3];
        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            dx[d] = dx[d] * inv_dr;
            dv[d] = home_v[d] - _mm256_set1_ps( other_v[d] );
        }

        auto rdotv = dx[0] * dv[0] + dx[1] * dv[1] + dx[2] * dv[2];

        auto sqrt_gammap = sqrt_diss_strength*wr;

        auto dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
        auto u = XorShift64Add( rng_state );
        for(int d=0; d<3; d++){
            assert(-0.5f <= u[d] && u[d] <= 0.5f );
        }
        auto randScale = sqrt_gammap * _mm256_set1_ps( scale_inv_sqrt_dt );
        auto randForce = randScale * u;

        auto scaled_force = _mm256_and_ps( (conForce+dissForce+randForce)*inv_dr, _mm256_castsi256_ps(active ));

        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            f[d] = scaled_force * dx[d];
        }

        return true;
    }

};

#endif