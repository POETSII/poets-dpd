#ifndef dpd_maths_primitives_hpp
#define dpd_maths_primitives_hpp

#include <cstdint>
#include <cstring>

inline double half(double x)
{ return x*(double)0.5; }

inline float half(float x)
{ return x*0.5f; }

inline double recip(double x)
{ return double(1.0)/x; }

inline float recip(float x)
{ return 1.0f/x; }

double recip_pow_half(double x);
float recip_pow_half(float x);

double pow_half(double x);
float pow_half(float x);

float absolute(float x);
double absolute(double x);

int floor_nn(double x);
int floor_nn(float x);

inline int round_impl(float x);

void memcpy32(uint32_t *a, const uint32_t *b, unsigned n);
void memzero32(uint32_t *a, const uint32_t *b, unsigned n);

template<unsigned N>
void memcpy32(uint32_t *a, const uint32_t *b);

template<unsigned N>
void memzero32(uint32_t *a);

// https://github.com/riscv-non-isa/riscv-toolchain-conventions
#if defined(__riscv) || defined (TINSEL)
#define PDPD_TINSEL
#include "dpd_maths_primitives_tinsel.hpp"
#else
#include "dpd_maths_primitives_native.hpp"
#endif

template<class T>
inline void memcpy32(T &a, const T &b)
{
    static_assert((sizeof(T)%4)==0);
    memcpy32<sizeof(T)/4>((uint32_t*)&a, (const uint32_t*)&b);
}

template<class T>
inline void memzero32(T &a)
{
    static_assert((sizeof(T)%4)==0);
    memzero32<sizeof(T)/4>((uint32_t*)&a);
}

#endif
