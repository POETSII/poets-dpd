#ifndef hash_hpp
#define hash_hpp

#include <cstdint>

inline uint64_t splitmix64(uint64_t z)
{
    z = (z ^ (z >> 30)) * uint64_t(0xBF58476D1CE4E5B9ull);
    z = (z ^ (z >> 27)) * uint64_t(0x94D049BB133111EBull);
    return z ^ (z >> 31);
}

uint32_t time_to_hash(double t, uint64_t seed)
{
    assert(t>=0);
    uint64_t tmp = (uint64_t)ldexp(t, 32);
    tmp += seed;
    return splitmix64(tmp) >> 32;
}

/*  This is a function which generates roughly random values
    for hash_rng(r,a,b), but with the constraint that hash_rng(t,a,b)==hash_rng(t,b,a)
    The return value should approximate 32 random bits, but the MSBs should be
    the priority.
*/
inline uint32_t hash_rng_sym(uint32_t t_hash, uint32_t a, uint32_t b)
{
    // TODO : This is terrible. Find the original version.
    uint32_t tmp=(t_hash^a)+(t_hash^b);
    tmp *= tmp>>16;
    tmp *= 0x1CE4E5B9ul;
    tmp ^= tmp>>16;
    tmp *= 0x1CE4E5B9ul;
    return tmp;
}

#endif