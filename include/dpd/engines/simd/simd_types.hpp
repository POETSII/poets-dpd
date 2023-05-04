#ifndef dpd_core_simd_types_hpp
#define dpd_core_simd_types_hpp

#if defined(__x86_64)
#include <immintrin.h>
#elif defined(__arm64)
#include <arm_neon.h>>
#else
#error "Unknown arch"
#endif

namespace simd_types
{

#ifdef __x86_64
    using simd_vecu32x8_t = __m256i;
    using simd_vecf32x8_t = __m256;

    static simd_vecu32x8_t simd_vecu32x8_make_zero()
    { return _mm256_setzero_si256(); }

    static simd_vecu32x8_t simd_vecu32x8_load(const uint32_t *src)
    { _mm256_loadu_si256( (const simd_vecu32x8_t *)src  ) }

#elif defined(__arm64)
    struct simd_vecu32x8_t{
        uint32x4_t lo, hi;
    };
    struct simd_vecf32x8_t{
        float32x4_t lo, hi;
    };

    inline simd_vecu32x8_t operator+(const simd_vecu32x8_t &a, const simd_vecu32x8_t &b)
    { return { vadd_u32(a.lo,a.hi), vadd_u32(a.hi,b.hi) }; }

    inline simd_vecf32x8_t operator+(const simd_vecf32x8_t &a, const simd_vecf32x8_t &b)
    { return { vadd_f32(a.lo,a.hi), vadd_f32(a.hi,b.hi) }; }

    inline simd_vecu32x8_t operator>(const simd_vecu32x8_t &a, const simd_vecu32x8_t &b)
    { return { vcgt_u32(a.lo,a.hi), vcgt_u32(a.hi,b.hi) }; }

    inline static simd_vecu32x8_t simd_vecu32x8_from_u64(
        uint64_t x3, uint64_t x2, uint64_t x1, uint64_t x0
    ){
        return {
            vreinterpretq_u32_u64(vcombine_u64(x0, x1)),
            vreinterpretq_u32_u64(vcombine_u64(x2,x3))
        };
    }

    inline  simd_vecu32x8_t simd_vecu32x8_from_u32(
        uint32_t x7, uint32_t x6, uint32_t x5, uint32_t x4,
        uint32_t x3, uint32_t x2, uint32_t x1, uint32_t x0
    ){
        return {
            vcombine_u32(
                vcreate_u32(x0|(uint64_t(x1)<<32)),
                vcreate_u32(x2|(uint64_t(x3)<<32))
            ),
            vcombine_u32(
                vcreate_u32(x4|(uint64_t(x5)<<32)),
                vcreate_u32(x6|(uint64_t(x7)<<32))
            )
        };
    }

    inline  simd_vecu32x8_t simd_vecu32x8_duplicate_scalar(uint32_t u)
    { return {vdupq_n_u32(u),vdupq_n_u32(u)}; }


    inline  simd_vecf32x8_t simd_vecf32x8_duplicate_scalar(float f)
    { return {vdupq_n_f32(f),vdupq_n_f32(f)}; }

    inline  simd_vecu32x8_t simd_vecu32x8_make_zero()
    { return {vdupq_n_u32(0),vdupq_n_u32(0)}; }

    inline  simd_vecf32x8_t simd_vecf32x8_make_zero()
    { return {vdupq_n_f32(0),vdupq_n_f32(0)}; }

    inline  simd_vecu32x8_t simd_vecu32x8_load(const uint32_t *src)
    { return { vld1_u32(src), vld1_u32(src+4) }; }

    inline  simd_vecf32x8_t simd_vecf32x8_load(const float *src)
    { return { vld1_f32(src), vld1_f32(src+4) }; }

#else
#error "Unknown architecture"
#endif

};

#endif
