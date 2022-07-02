#include "dpd/core/non_conflict_grid.hpp"

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>

#include <time.h>

#include <iostream>

__attribute__((noinline))
int get_reps(int x, int y, int z)
{
    return 100000 + (x&1)*10000 + (y&1)*10000 + (z&1)*10000;
};

__attribute__((noinline))
int get_inc()
{
    return 1;
}

template<class TImpl>
void test(int n)
{
    std::vector<int> counters(n*n*n, 0);

    auto make_task=[&](int x, int y, int z)
    {
        return std::function<void()>([&,x,y,z](){
            volatile int *c = &counters[x+y*n+z*n*n];
            int acc=*c;
            for(int i=0; i<get_reps(x,y,z); i++){
                acc=acc+get_inc();
            }
            *c=acc;
        });
    };

    TImpl impl(n,n,n,make_task);

    for(int i=0; i<100; i++){
        impl.execute();
    }

    for(int i=0; i<n*n*n; i++){
        int e=get_reps( i%n, (i/n)%n, i/(n*n) ) * 100;
        if(counters[i]!=e){
            std::cerr<<"Broken : e="<<e<<", got="<<counters[i]<<"\n";
            exit(1);
        }
    }
}

double now()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

int main()
{
    int max_parallelism=tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    if(getenv("PDPD_NUM_THREADS")){
	    max_parallelism=std::atoi(getenv("PDPD_NUM_THREADS"));
    }
    std::cerr<<"TBB is using "<<max_parallelism<<" threads.\n";
    tbb::global_control tbb_control_threads(tbb::global_control::max_allowed_parallelism, max_parallelism);

    for(int i=2; i<32; i+=2){
        std::cerr<<"i="<<i<<"\n";
        double a=now();
        test<graph_grid_by_flow>(i);
        double b=now();
        test<graph_grid_by_parfor>(i);
        double c=now();
        std::cerr<<"i="<<i<<", flow="<<b-a<<", parfor="<<c-b<<"\n";
    }
}
