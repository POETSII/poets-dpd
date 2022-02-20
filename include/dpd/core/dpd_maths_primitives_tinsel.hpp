#ifndef dpd_maths_primitives_tinsel_hpp
#define dpd_maths_primitives_tinsel_hpp

#include <cstdint>
#include <cstring>

// inverse square root using bit twiddling cleverness (taken from http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf)
// uses newton raphson
inline float recip_pow_half(float x) {
    float xhalf = 0.5f*x;

    union{
        float f;
        int32_t i;
        uint32_t u;
    } tmp;
    tmp.f=x;
    tmp.i = 0x5f3759df - (tmp.i>>1); // gives initial guess y0
    x = tmp.f; // convert bits back to float
    // Relative error is:
    // - 0.175228 after 1 iteration
    // - 4.66e-004 after 2 iterations
    // - (estimated) 2^-22 after 3 iterations
    for(int i=0; i<3; i++) {
        x = x*(1.5f-xhalf*x*x); // Newton step, repeating increases accuracy
    }
    return x;
}

// sqrt using the inverse square root function
inline float pow_half(float x) {
   //return 1.0/inv_sqrt(x);

   // https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
   // We have more integer than float isue bandwidth, plus fast divides (well, same speed as anything else),
   // so might as well use this one.

   //tinselAssert(x>=0);

   union
   {
       int i;
       float x;
   } u;
   u.x = x;
   u.i = (1<<29) + (u.i >> 1) - (1<<22);

   // Two Babylonian Steps (simplified from:)
   // u.x = 0.5f * (u.x + x/u.x);
   // u.x = 0.5f * (u.x + x/u.x);

   //u.x =       u.x + x/u.x;
   //u.x = 0.25f*u.x + x/u.x;

   u.x = 0.5f * (u.x + x/u.x);
   u.x = 0.5f * (u.x + x/u.x);
   u.x = 0.5f * (u.x + x/u.x);

#ifdef DOUBLE_SQRT
   u.x = 0.5f * (u.x + x/u.x);
   u.x = 0.5f * (u.x + x/u.x);
   u.x = 0.5f * (u.x + x/u.x);
#endif

   return u.x;
}

inline float absolute(float x)
{
    union{
        float f;
        int32_t i;
        uint32_t u;
    } tmp;
    tmp.f=x;
    tmp.u &= 0x7FFFFFFFul;
    return tmp.f;
}

inline int floor_nn(float x)
{
    // WARNING : The round-to-even behaviour of tinsel is a problem here
    // I think this is equivalent to floor, _assumign_ that int conversion is round to even
    const float delta=0.499999970198f;
    int r=(int)(x-delta); // Warning! This will get converted to round to nearest even
    
    //assert( r= floorf(x) ); // Note: this is _expected to fail in x86. It should pass in tinsel though.

    return r;
}

void memcpy32(uint32_t *a, const uint32_t *b, unsigned n);
void memzero32(uint32_t *a, unsigned n);

#define TINSEL_MEMCPY_CROSSOVER 4

template<unsigned N>
inline void memcpy32(uint32_t *a, const uint32_t *b)
{
    if constexpr (N<=TINSEL_MEMCPY_CROSSOVER){
        #pragma GCC unroll 4
        for(unsigned i=0; i<N; i++){
            a[i] = b[i];
        }
    }else{
        memcpy32(a,b,N);
    }
}

#define TINSEL_MEMZERO_CROSSOVER 8

template<unsigned N>
void memzero32(uint32_t *a)
{
    if constexpr (N<=TINSEL_MEMZERO_CROSSOVER){
        #pragma GCC unroll 8
        for(unsigned i=0; i<N; i++){
            a[i] = 0;
        }
    }else{
        memzero32(a,N);
    }
}


inline int round_impl(float x)
{
    return (int)x; // Tinsel only does round nearest even
}


#endif