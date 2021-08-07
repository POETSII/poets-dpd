#include "basic_dpd_engine_v5_raw_tinsel.hpp"

#include <tinsel.h>
#include <POLite.h>

using Thread = typename BasicDPDEngineV5RawTinsel<POLiteHW<>>::Thread;

float absolute(float x)
{
    uint32_t y=*(uint32_t*)&x;
    y &= 0x7FFFFFFFul;
    return *(float*)&y;
}

// inverse square root using bit twiddling cleverness (taken from http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf)
// uses newton raphson
float recip_pow_half(float x) {
    float xhalf = 0.5f*x;
    int i = *(int*)&x; // get bits for floating value
    i = 0x5f3759df - (i>>1); // gives initial guess y0
    x = *(float*)&i; // convert bits back to float
    for(int i=0; i<5; i++) {
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


void * memcpy ( void * destination, const void * source, size_t num )
{
    for(size_t i=0; i<num; i++){
        ((char*)destination)[i] = ((const char*)source)[i];
    }
    return destination;
}

inline int floor_nn(float x)
{
    // MASSIVE TODO : The round-to-even behaviour of tinsel is a problem here
    int r=0;
    while(x>=1.0f){
        x=x-1;
        r=r+1;
    }
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
