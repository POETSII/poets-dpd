#include "basic_dpd_engine_v4_raw_handlers.hpp"

using Handlers=BasicDPDEnginev4RawHandlers<>;

int main()
{}

extern "C" void *memcpy(void *__restrict dst, const void *__restrict src, size_t n)
{
    for(size_t i=0; i<n; i++){
        ((char*)dst)[i] = ((const char*)src)[i];
    }
    return dst;
}


__attribute__((externally_visible)) void *x0=(void*)Handlers::on_barrier;

__attribute__((externally_visible)) void *x1=(void*)Handlers::on_send_migrate;
__attribute__((externally_visible)) void *x2=(void*)Handlers::on_recv_migrate;


__attribute__((externally_visible)) void *x3=(void*)Handlers::on_send_force;
__attribute__((externally_visible)) void *x4=(void*)Handlers::on_recv_force;
__attribute__((externally_visible)) void *x5=(void*)Handlers::on_send_share;
__attribute__((externally_visible)) void *x6=(void*) Handlers::on_recv_share;