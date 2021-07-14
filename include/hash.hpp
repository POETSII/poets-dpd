#ifndef hash_hpp
#define hash_hpp

#include <cstdint>

    inline uint64_t splitmix64(uint64_t z)
    {
        z = (z ^ (z >> 30)) * uint64_t(0xBF58476D1CE4E5B9ull);
        z = (z ^ (z >> 27)) * uint64_t(0x94D049BB133111EBull);
        return z ^ (z >> 31);
    }

#endif