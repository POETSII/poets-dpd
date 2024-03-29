#ifndef hash_hpp
#define hash_hpp

#include <cstdint>
#include <algorithm>
#include <cassert>

#if __x86_64__ 
#include <immintrin.h>
#endif

inline uint64_t splitmix64(uint64_t z)
{

    z = (z ^ (z >> 30)) * uint64_t(0xBF58476D1CE4E5B9ull);
    z = (z ^ (z >> 27)) * uint64_t(0x94D049BB133111EBull);
    return z ^ (z >> 31);
}

inline uint64_t riscv_mix64_m2(uint64_t x)
{
    // From custom hash prospector. Not as good as splitmix64 in quality, but much faster/smaller on riscv
    const uint32_t C=0xd392d2a7; // 2
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    return x;
}

#if 0
// Does hash on 4 64-bit integers
inline __m256i riscv_mix64_m2_v4q(__m256i x)
{
    const __m256i C=_mm256_set1_epi64x (0xd392d2a7); // 2
    x = x ^ _mm256_srli_epi64(x,32);     // 1
    x = x * C;           // 4
    x = x ^ _mm256_srli_epi64(x,32);     // 1
    x = x * C;           // 4
    x = x ^ _mm256_srli_epi64(x,32);     // 1
    return x;
}
#endif

inline uint64_t riscv_mix64_m3(uint64_t x)
{
    // From custom hash prospector. As good as splitmix64, but faster/smaller on riscv
    const uint32_t C=0xed85aebf; // 2
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    return x;
}


inline uint64_t get_t_hash(uint32_t t, uint64_t seed)
{
    uint64_t base=seed + riscv_mix64_m2(t);
    return riscv_mix64_m2(base);
}

/*  This is a function which generates roughly random values
    for hash_rng(r,a,b), but with the constraint that hash_rng(t,a,b)==hash_rng(t,b,a)
    The return value should approximate 32 random bits, but the MSBs should be
    the priority.
*/
inline uint32_t hash_rng_sym_good(uint64_t t_hash, uint32_t a, uint32_t b)
{
    if(a>b){
        std::swap(a,b);
    }
    // m3 is safer from randomness testing perspective, but m2 is ok in practise for many beads.
    return riscv_mix64_m3( (a|(uint64_t(b)<<32)) ^ t_hash);
}

inline uint32_t hash_rng_sym(uint64_t t_hash, uint32_t a, uint32_t b)
{
    if(a>b){
        std::swap(a,b);
    }
    // m3 is safer from randomness testing perspective, but m2 is ok in practise for many beads.
    return riscv_mix64_m2( (a|(uint64_t(b)<<32)) ^ t_hash); // Return LSBS
}

inline uint32_t hash_rng_sym_crappy(uint64_t t_hash, uint32_t a, uint32_t b)
{
    assert(t_hash & 0x1);
    assert(t_hash & 0x100000000ull);

    uint32_t aa=0xed85aebfull*a + 0xed85aebfull*b + (t_hash&0xFFFFFFFFul);
    aa ^= aa>>16;
    aa += a+b+(t_hash>>32); 
    return aa * 0xed85aebful;
}

// Calculate checksum except for the last 4 bytes (assumed to be the checksum field).
    template<class T>
    static uint32_t calc_checksum(const T &x)
    {
        static_assert( (sizeof(T) % 4) == 0 );
        const uint32_t *src=(const uint32_t*)&x;
        uint32_t res=0;
        for(unsigned i=0; i<sizeof(T)/4-1; i++){
            res = res + (res>>16) + src[i];
        }
        return res;
    }

#endif