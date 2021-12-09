#include "dpd/engines/basic/basic_dpd_engine_v7_raw_tinsel.hpp"

#include <tinsel.h>
#include <POLite.h>

#include "POLiteHW.h"

constexpr bool USE_X_CACHE=true;
using Thread = typename BasicDPDEngineV7RawTinsel<POLiteHW<>,USE_X_CACHE>::Thread;

using f2i_t = union{
        float f;
        int32_t i;
        uint32_t u;
    };

float absolute(float x)
{
    f2i_t tmp;
    tmp.f=x;
    tmp.u &= 0x7FFFFFFFul;
    return tmp.f;
}

// inverse square root using bit twiddling cleverness (taken from http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf)
// uses newton raphson
float recip_pow_half(float x) {
    float xhalf = 0.5f*x;

    f2i_t tmp;
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
float pow_half(float x) {
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


inline int floor_nn(float x)
{
    // WARNING : The round-to-even behaviour of tinsel is a problem here
    // I think this is equivalent to floor, _assumign_ that int conversion is round to even
    const float delta=0.499999970198f;
    int r=(int)(x-delta); // Warning! This will get converted to round to nearest even
    
    assert( r= floorf(x) ); // Note: this is _expected to fail in x86. It should pass in tinsel though.

    return r;
}

int main()
{
  // Point thread structure at base of thread's heap
  Thread* thread = (Thread*) tinselHeapBaseSRAM();
  
  // Invoke interpreter
  thread->run();

  return 0;
}
