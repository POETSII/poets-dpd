#ifndef dpd_maths_core_simd_hpp
#define dpd_maths_core_simd_hpp

#include <immintrin.h>

#include "dpd/core/hash.hpp"
#include "dpd/maths/dpd_maths_core.hpp"
#include "dpd/core/logging.hpp"

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
    static __m256 XorShift64Add(__m256i &x)
    {
        __m256i a=x;

        x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
        x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
        x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));

        __m256i b=_mm256_shuffle_epi32(x, (2<<6)|(3<<4)|(0<<2)|(1<<0));

        auto u= _mm256_add_epi32(a,b);

        const __m256 scale = _mm256_set1_ps(0.00000000023283064365386962890625f);
        __m256 uf = _mm256_cvtepi32_ps(u);
        return uf*scale;
    }

    static __m256 POETSHashV1(uint64_t base, __m256i id1, __m256i id2)
    {
        uint32_t id1b[8];
        uint32_t id2b[8];

        _mm256_storeu_si256((__m256i*)id1b,id1);
        _mm256_storeu_si256((__m256i*)id2b, id2);

        uint32_t hash[8];
        for(unsigned i=0; i<8; i++){
            hash[i]=hash_rng_sym(base, id1b[i], id2b[i]);
        }

        auto u= _mm256_loadu_si256((__m256i*)hash);

        const __m256 scale = _mm256_set1_ps(0.00000000023283064365386962890625f);
        __m256 uf = _mm256_cvtepi32_ps(u);
        return uf*scale;
    }


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

    // Add together the 4 upper and 4 lower and return two sums
    static inline std::pair<float,float> mm256_reduce_add_ps_dual_vec4(__m256 x) {
        __m256 s0 = _mm256_shuffle_ps(x, x, (3<<4) | (1<<0) );    // ( -, x7, -, x5, -, x3, -, x1)
        __m256 s1 = _mm256_add_ps( x, s0 );                       // (-, x7+x6, -, x5+x4, -, x3+x2, -, x1+x0)

        __m256 s2 = _mm256_shuffle_ps(s1, s1, (2<<0) );    // ( -, -, -, x7+x6, -, -, -, x3+x2)
        __m256 s3 = _mm256_add_ps(s1, s2);               // (-, -, -, x7+x6+x5+x4, -, -, - x3+x2+x1+x0)

        float lo=s3[0];
        float hi=s3[4];

        return {lo,hi};
    }

    enum Flags
    {
        Flag_EnableLogging  =1,
        
        Flag_RngHash           =2,
        Flag_RngXorShift64Add  =4,
        Flag_RngZero  =8,

        Flag_RngMask           =2|4|8
    };

    /* Returns true if there is any interaction (non-zero force), false otherwise.
        The f values are only valid iff it returns true. f is not set to zero
        if there are no interactions.

        This function treats all home lanes as active. It is up to the caller to
        make sure they are valid (e.g. by moving inactive lanes to an x value which
        can't interact) or by applying masking externally.

        ACTIVE_BEADS_FOR_LOGGING can be used to control how many beads are considered
        active for logging purposes. It does not change the calculations.
    */
    template<Flags TFlags, unsigned MAX_BEAD_TYPES, unsigned ACTIVE_BEADS_FOR_LOGGING=8>
    static bool interact_vec8_to_scalar(
        const float scale_inv_sqrt_dt,
        const float conservative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        const float sqrt_dissipative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        
        uint64_t t_hash,
        __m256i &rng_state,

        const __m256i home_bead_ids, 
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],

        uint32_t other_bead_id,
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

        __m256 dr, inv_dr;
        if(0){
            dr = _mm256_sqrt_ps (dr2);
            inv_dr = 1.0f / dr;
        }else if(0){
            inv_dr = _mm256_rsqrt_ps(dr2);
            dr = inv_dr * dr2;
        }else{
            // This seems to be marginally faster. They both go down the
            // same unit, but there are no dependencies
            inv_dr = _mm256_rsqrt_ps(dr2);
            dr = _mm256_sqrt_ps(dr2);
        }

        __m256 con_strength_row = _mm256_loadu_ps( conservative_matrix + 8 * other_bead_type);
        __m256 sqrt_diss_strength_row = _mm256_loadu_ps( sqrt_dissipative_matrix + 8 * other_bead_type);
        __m256 con_strength = _mm256_permutevar8x32_ps( con_strength_row, home_bead_types );
        __m256 sqrt_diss_strength = _mm256_permutevar8x32_ps( sqrt_diss_strength_row, home_bead_types );

        auto wr = 1.0f - dr;
        auto conForce = con_strength * wr;

        __m256 dv[3];
        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            dv[d] = home_v[d] - _mm256_set1_ps( other_v[d] );
        }

        auto rdotv = dx[0] * dv[0] + dx[1] * dv[1] + dx[2] * dv[2];

        auto sqrt_gammap = sqrt_diss_strength*wr;

        auto dissForce = -sqrt_gammap*sqrt_gammap*rdotv*inv_dr;
        __m256 u;
        if(TFlags & Flag_RngXorShift64Add){
            u=XorShift64Add( rng_state );
        }else if (TFlags & Flag_RngHash){
            u=POETSHashV1(t_hash, home_bead_ids, _mm256_set1_epi32(other_bead_id));
        }else if (TFlags & Flag_RngZero){
            u=_mm256_setzero_ps();
        }else{
            throw std::logic_error("No Rng method selected at compile-time.");
        }
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

        if(TFlags & Flag_EnableLogging && ForceLogging::logger()){
            BeadHash thash{other_bead_id};
            uint32_t activeb[8], home_bead_idsb[8];
            _mm256_storeu_si256((__m256i*)activeb, active);
            _mm256_storeu_si256((__m256i*)home_bead_idsb, home_bead_ids);
            for(unsigned i=0; i<ACTIVE_BEADS_FOR_LOGGING; i++){
                if(!activeb[i]){
                    continue;
                }
                BeadHash hhash{home_bead_idsb[i]};
                double ddx[3]={dx[0][i],dx[1][i],dx[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dx", 3,ddx);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dr", 1,&dr[i]);

                double ddv[3]={dv[0][i],dv[1][i],dv[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dv", 3,ddv);
                double dd=sqrt_diss_strength[i] * sqrt_diss_strength[i];
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-diss-strength", 1, &dd);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-invrootdt", 1, &scale_inv_sqrt_dt);
                double gammap=sqrt_gammap[i]*sqrt_gammap[i];
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-gammap", 1, &gammap);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rng", 1, &u[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-con", 1, &conForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-diss", 1,&dissForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rng-scale",1, &randScale[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rand",1, &randForce[i]);
                double ff[3]={f[0][i],f[1][i],f[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"f_next_dpd", 3,ff);

                for(int d=0; d<3; d++){ ddx[d]=-ddx[d]; ddv[d]=-ddv[d]; }
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dx", 3,ddx);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dr", 1,&dr[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dv", 3,ddv);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-diss-strength", 1, &dd);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-invrootdt", 1, &scale_inv_sqrt_dt);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-gammap", 1, &gammap);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rng", 1, &u[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-con", 1, &conForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-diss", 1,&dissForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rng-scale",1, &randScale[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rand",1, &randForce[i]);
                for(int d=0; d<3; d++){ ff[d] = -ff[d]; }
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"f_next_dpd", 3,ff);
            }
        }

        return true;
    }

    template<Flags TFlags, unsigned MAX_BEAD_TYPES>
    static bool interact_dual_vec4_to_dual_scalar(
        const float scale_inv_sqrt_dt,
        const float conservative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        const float sqrt_dissipative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        
        uint64_t t_hash,
        __m256i &rng_state,
        
        const __m256i home_bead_ids,
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],

        uint32_t otherA_bead_id,
        unsigned otherA_bead_type,
        float otherA_x[3],
        float otherA_v[3],

        uint32_t otherB_bead_id,
        unsigned otherB_bead_type,
        float otherB_x[3],
        float otherB_v[3],

        __m256 f[3] // Force on home beads
    ){
        static_assert(MAX_BEAD_TYPES <= 8);

        __m256 dx[3], ox[3];
        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            ox[d]=_mm256_set1_ps(otherA_x[d]);
            ox[d]=_mm256_insertf128_ps(ox[d], _mm_set1_ps(otherB_x[d]), 1);

            dx[d] = home_x[d] - ox[d];
        }
        __m256 dr2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        __m256i active = (dr2 > 0.00001f && 1.0f > dr2);
        if(0==_mm256_movemask_epi8(active)){
           return false;
        }

        __m256 dr, inv_dr;
        if(0){
            dr = _mm256_sqrt_ps (dr2);
            inv_dr = 1.0f / dr;
        }else if(0){
            inv_dr = _mm256_rsqrt_ps(dr2);
            dr = inv_dr * dr2;
        }else{
            // This seems to be marginally faster. They both go down the
            // same unit, but there are no dependencies
            inv_dr = _mm256_rsqrt_ps(dr2);
            dr = _mm256_sqrt_ps(dr2);
        }

        __m256i other_bead_type = _mm256_set1_epi32(otherA_bead_type);
        other_bead_type = _mm256_insertf128_si256(other_bead_type,_mm_set1_epi32(otherB_bead_type),1);

        __m256 con_strength_row = _mm256_loadu_ps( conservative_matrix + 8 * otherA_bead_type);
        __m256 sqrt_diss_strength_row = _mm256_loadu_ps( sqrt_dissipative_matrix + 8 * otherA_bead_type);

        // Not sure if the branch cost is worth it
        if(otherA_bead_type!=otherB_bead_type){
            __m128 con_strength_B_row = _mm_loadu_ps( conservative_matrix + 8 * otherB_bead_type );
            __m128 sqrt_diss_strength_B_row = _mm_loadu_ps( sqrt_dissipative_matrix + 8 * otherB_bead_type );

            con_strength_row = _mm256_insertf128_ps(con_strength_row, con_strength_B_row, 1);
            sqrt_diss_strength_row = _mm256_insertf128_ps(con_strength_row, con_strength_B_row, 1);
        }

        __m256 con_strength = _mm256_permutevar8x32_ps( con_strength_row, home_bead_types );
        __m256 sqrt_diss_strength = _mm256_permutevar8x32_ps( sqrt_diss_strength_row, home_bead_types );
        
        auto wr = 1.0f - dr;
        auto conForce = con_strength * wr;

        __m256 dv[3];
        #pragma GCC unroll(3)
        for(int d=0; d<3; d++){
            __m256 ov=_mm256_set1_ps(otherA_v[d]);
            ov=_mm256_insertf128_ps(ov, _mm_set1_ps(otherB_v[d]), 1);

            dv[d] = home_v[d] - ov;
        }

        auto rdotv = dx[0] * dv[0] + dx[1] * dv[1] + dx[2] * dv[2];

        auto sqrt_gammap = sqrt_diss_strength*wr;

        auto dissForce = -sqrt_gammap*sqrt_gammap*rdotv*inv_dr;
        __m256 u;
        if(TFlags & Flag_RngXorShift64Add){
            u=XorShift64Add( rng_state );
        }else if (TFlags & Flag_RngHash){
            __m256i other_bead_id=_mm256_insertf128_si256(_mm256_set1_epi32(otherA_bead_id), _mm_set1_epi32(otherB_bead_id), 1);
            u=POETSHashV1(t_hash, home_bead_ids, other_bead_id);
        }else if (TFlags & Flag_RngZero){
            u=_mm256_setzero_ps();
        }else{
            throw std::logic_error("No Rng method selected at compile-time.");
        }
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

        if(TFlags & Flag_EnableLogging && ForceLogging::logger()){
            uint32_t activeb[8], home_bead_idsb[8];
            _mm256_storeu_si256((__m256i*)activeb, active);
            _mm256_storeu_si256((__m256i*)home_bead_idsb, home_bead_ids);


            for(unsigned i=0; i<8; i++){
                if(!activeb[i]){
                    continue;
                }
                BeadHash hhash{home_bead_idsb[i]};
                BeadHash thash{ i>=4 ? otherB_bead_id : otherA_bead_id };
                assert(hhash!=thash);
                double ddx[3]={dx[0][i],dx[1][i],dx[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dx", 3,ddx);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dr", 1,&dr[i]);
                double ddv[3]={dv[0][i],dv[1][i],dv[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dv", 3,ddv);
                double dd=sqrt_diss_strength[i] * sqrt_diss_strength[i];
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-diss-strength", 1, &dd);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-invrootdt", 1, &scale_inv_sqrt_dt);
                double gammap=sqrt_gammap[i]*sqrt_gammap[i];
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-gammap", 1, &gammap);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rng", 1, &u[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-con", 1, &conForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-diss", 1,&dissForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rng-scale",1, &randScale[i]);
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"dpd-rand",1, &randForce[i]);
                double ff[3]={f[0][i],f[1][i],f[2][i]};
                ForceLogging::logger()->LogBeadPairProperty(hhash,thash,"f_next_dpd", 3,ff);

                for(int d=0; d<3; d++){ ddx[d]=-ddx[d]; ddv[d]=-ddv[d]; }
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dx", 3,ddx);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dr", 1,&dr[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dv", 3,ddv);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-diss-strength", 1, &dd);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-invrootdt", 1, &scale_inv_sqrt_dt);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-gammap", 1, &gammap);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rng", 1, &u[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-con", 1, &conForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-diss", 1,&dissForce[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rng-scale",1, &randScale[i]);
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"dpd-rand",1, &randForce[i]);
                for(int d=0; d<3; d++){ ff[d] = -ff[d]; }
                ForceLogging::logger()->LogBeadPairProperty(thash,hhash,"f_next_dpd", 3,ff);
            }
        }

        return true;
    }
};

#endif