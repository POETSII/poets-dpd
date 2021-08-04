#ifndef dpd_maths_primitives_hpp
#define dpd_maths_primitives_hpp

#include <cstdint>

inline double half(double x)
{ return x*0.5; }

inline float half(float x)
{ return x*0.5f; }

inline double recip(double x)
{ return 1.0/x; }

inline float recip(float x)
{ return 1.0f/x; }

// These wierd names are to stop magic from the C library and compiler getting in the way

inline double recip_pow_half(double x)
#ifdef TINSEL
;
#else
{ return 1.0/sqrt(x); }
#endif

inline float recip_pow_half(float x)
#ifdef TINSEL
;
#else
{ return 1.0f/sqrtf(x); }
#endif

inline double pow_half(double x)
#ifdef TINSEL
;
#else
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




#endif