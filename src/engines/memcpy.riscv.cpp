
#include <cstdint>

extern "C"  void * memcpy ( void * destination, const void * source, uint32_t num ) 
{
    for(uint32_t i=0; i<num; i++){
        ((char*)destination)[i] = ((const char*)source)[i];
    }
    return destination;
}

void /*__attribute__((noinline)) __attribute__((optimize("no-tree-loop-distribute-patterns")))*/ memcpy32(uint32_t *a, const uint32_t *b, unsigned n)
{
    for(unsigned i=0; i<n; i++){
        a[i]=b[i];
    }
}

void /*__attribute__((noinline)) __attribute__((optimize("no-tree-loop-distribute-patterns")))*/ memzero32(uint32_t *a, unsigned n)
{
    for(unsigned i=0; i<n; i++){
        a[i]=0;
    }
}
