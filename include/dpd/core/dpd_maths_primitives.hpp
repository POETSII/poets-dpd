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

// These wierd names are to stop magic from the C library and compiler getting in the way

inline double recip_pow_half(double x)
#ifdef TINSEL
;
#else
{ return 1.0/sqrt(x); }
#endif


#ifdef TINSEL
float recip_pow_half(float x);
#else
inline float recip_pow_half(float x)
{ return 1.0f/sqrtf(x); }
#endif


#ifdef TINSEL
double pow_half(double x);
#else
inline double pow_half(double x)
{ return sqrt(x); }
#endif

inline float pow_half(float x)
#ifdef TINSEL
;
#else
{ return sqrtf(x); }
#endif

inline float absolute(float x)
#ifdef TINSEL
;
#else
{ return fabsf(x); }
#endif

inline double absolute(double x)
#ifdef TINSEL
;
#else
{ return fabs(x); }
#endif

inline int floor_nn(double x)
#ifdef TINSEL
;
#else
{ return (int)floor(x); }
#endif

inline int floor_nn(float x)
#ifdef TINSEL
;
#else
{ return (int)floor(x); }
#endif

#ifdef TINSEL
void memcpy32(uint32_t *a, const uint32_t *b, unsigned n);
#else
inline void memcpy32(uint32_t *a, const uint32_t *b, unsigned n)
{
    memcpy(a, b, n*4);
}
#endif


#endif