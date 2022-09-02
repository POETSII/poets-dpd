#ifndef dpd_maths_primitives_native_hpp
#define dpd_maths_primitives_native_hpp

#include <cstdint>
#include <cstring>
#include <cmath>

inline double recip_pow_half(double x)
{ return 1.0/sqrt(x); }

inline float recip_pow_half(float x)
{ return 1.0f/sqrtf(x); }


inline double pow_half(double x)
{ return sqrt(x); }

inline float pow_half(float x)
{ return sqrtf(x); }


inline float absolute(float x)
{ return fabsf(x); }

inline double absolute(double x)
{ return fabs(x); }


inline int floor_nn(double x)
{ return (int)floor(x); }

inline int floor_nn(float x)
{ return (int)floor(x); }


inline void memcpy32(uint32_t *a, const uint32_t *b, unsigned n)
{
    memcpy(a, b, n*4);
}

inline void memcpy32(uint32_t *a, const volatile uint32_t *b, unsigned n)
{
    memcpy32(a, (const uint32_t*)b, n);
}

inline void memcpy32(volatile uint32_t *a, const uint32_t *b, unsigned n)
{
    memcpy32((uint32_t*)a, (const uint32_t*)b, n);
}

inline void memzero32(uint32_t *a, unsigned n)
{
    memset(a, 0, n*4);
}

template<unsigned N>
void memcpy32(uint32_t *a, const uint32_t *b)
{
    memcpy32(a,b,N);
}

template<unsigned N>
void memcpy32(uint32_t *a, const volatile uint32_t *b)
{
    memcpy32(a,b,N);
}

template<unsigned N>
void memcpy32(volatile uint32_t *a, const uint32_t *b)
{
    memcpy32(a,b,N);
}

template<unsigned N>
void memzero32(uint32_t *a)
{
    memzero32(a,N);
}


inline void memswap32(uint32_t *a, uint32_t *b, unsigned n)
{
    std::swap_ranges(a, a+n, b);
}

template<unsigned N>
inline void memswap32(uint32_t *a, uint32_t *b)
{
    memswap32(a, b, N);
}

inline int round_impl(float x)
{
    return roundf(x);
}

inline uint32_t tinselCycleCount()
{
    return 0;
}

#endif