#include "dpd/core/non_conflict_grid.hpp"

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>

#include <time.h>

#include <iostream>
#include <cmath>

double now()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

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

struct timing_info_t
{
    double create;
    double first;
    double mean;
    double stddev;
};

std::ostream &operator<<(std::ostream &dst, const timing_info_t &ti)
{   
    return dst<<"init="<<ti.create<<", first="<<ti.first<<", mean="<<ti.mean<<", stddev="<<ti.stddev;
}

template<class TImpl>
timing_info_t test(int n)
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

    double t1=now();
    TImpl impl(n,n,n,make_task);

    int reps=100;
    
    double t2=now();
    impl.execute();
    double t3=now();
    double sumSqr = std::pow(t3-t2,2);

    double tStart=t3;
    for(int i=1; i<reps; i++){
        impl.execute();
        double tEnd=now();
        sumSqr += std::pow(tEnd-tStart,2);
        tStart=tEnd;
    }
    double t4=now();

    for(int i=0; i<n*n*n; i++){
        int e=get_reps( i%n, (i/n)%n, i/(n*n) ) * 100;
        if(counters[i]!=e){
            std::cerr<<"Broken : e="<<e<<", got="<<counters[i]<<"\n";
            exit(1);
        }
    }

    double mean=(t4-t2)/reps;
    double stddev=std::sqrt(sumSqr/reps - mean*mean);
    return {t2-t1, t3-t2, mean, stddev};
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
        std::cerr<<"flow, "<<test<graph_grid_by_flow>(i)<<"\n";
        std::cerr<<"parfor, "<<test<graph_grid_by_parfor>(i)<<"\n";
        std::cerr<<"affinity, "<<test<graph_grid_by_affinity_queue>(i)<<"\n";
        
    }
}
