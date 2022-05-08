
#include <cstdint>

extern "C" void __attribute__((noinline)) __attribute__((optimize("no-tree-loop-distribute-patterns"))) memcpy32(volatile uint32_t *a, const volatile uint32_t *b, unsigned n)
{
    for(unsigned i=0; i<n; i++){
        a[i]=b[i];
    }
}

extern "C"  void * __attribute__((noinline)) __attribute__((optimize("no-tree-loop-distribute-patterns"))) memcpy ( void * destination, const void * source, uint32_t num ) 
{
    volatile char *dst=(char *)destination;
    volatile const char *src=(const char *)source;

    for(uint32_t i=0; i<num; i++){
        dst[i]=src[i];
    }
    return destination;
}

void /*__attribute__((noinline)) __attribute__((optimize("no-tree-loop-distribute-patterns")))*/ memzero32(uint32_t *a, unsigned n)
{
    for(unsigned i=0; i<n; i++){
        a[i]=0;
    }
}
