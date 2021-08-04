#ifndef hash_hpp
#define hash_hpp

#include <cstdint>
#include <algorithm>
#include <cassert>

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

/*
a=a*C;
b=b*C;
a ^= (b>>16);
b ^= (a>>16);
a=a*C;
b=b*C;
a ^= (b>>16);
b ^= (a>>16);
*/



inline uint64_t next_t_hash(uint64_t &seed)
{
    seed += 0x9E3779B97F4A7C15ull;
    return riscv_mix64_m2(seed) | 0x0000000100000001ull;
}

/*  This is a function which generates roughly random values
    for hash_rng(r,a,b), but with the constraint that hash_rng(t,a,b)==hash_rng(t,b,a)
    The return value should approximate 32 random bits, but the MSBs should be
    the priority.
*/
inline uint32_t hash_rng_sym(uint64_t t_hash, uint32_t a, uint32_t b)
{
    // TODO : This is terrible. Find the original version.
    if(a>b){
        std::swap(a,b);
    }
    return riscv_mix64_m3( (a|(uint64_t(b)<<32)) ^ t_hash);
}

inline uint32_t hash_rng_sym_new(uint64_t t_hash, uint32_t a, uint32_t b)
{
    assert(t_hash & 0x1);
    assert(t_hash & 0x100000000ull);

    auto la = uint32_t(t_hash&0xFFFFFFFFul) * (a+b);
    auto lb = uint32_t(t_hash>>32) * (a^b);
    uint32_t tmp = la^lb;
    return tmp;
}

#endif