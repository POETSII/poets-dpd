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


uint64_t next_t_hash(uint64_t &seed)
{
    seed += 0x9E3779B97F4A7C15ull;
    return splitmix64(seed) | 0x0000000100000001ull;
}

/*  This is a function which generates roughly random values
    for hash_rng(r,a,b), but with the constraint that hash_rng(t,a,b)==hash_rng(t,b,a)
    The return value should approximate 32 random bits, but the MSBs should be
    the priority.
*/
inline uint32_t hash_rng_sym_old(uint32_t t_hash, uint32_t a, uint32_t b)
{
    // TODO : This is terrible. Find the original version.
    uint32_t la=std::min(a,b);
    uint32_t lb=std::max(a,b);
    uint64_t lla=splitmix64(la|(uint64_t(lb)<<32));
    uint64_t llb=splitmix64(lla^(t_hash^(uint64_t(t_hash)<<32)));
    return llb;
}

inline uint32_t hash_rng_sym(uint64_t t_hash, uint32_t a, uint32_t b)
{
    assert(t_hash & 0x1);
    assert(t_hash & 0x100000000ull);

    auto la = (t_hash&0xFFFFFFFFul) * (a+b);
    auto lb = (t_hash>>32) * (a^b);
    uint32_t tmp = la^lb;
    return tmp;
}

#endif