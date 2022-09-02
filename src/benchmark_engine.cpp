#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"
#include <tbb/global_control.h>

#include <random>
#include <iostream>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>


void usage()
{
    fprintf(stderr, "benchmark_engine : engine-name mode \n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr, "  modes:\n");
    fprintf(stderr, "    uniform-n : Box full of beads with density 3 and volume n^3 \n");
    exit(1);
}

double now()
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

WorldState make_uniform(int dim)
{
    WorldStateBuilder b({dim,dim,dim});
    WorldState &s=b.data();

    s.dt=0.01;
    s.interactions.push_back({1,0.1});
    s.lambda=0.5;
    s.t=0;
    b.add_bead_type("W");
    b.add_polymer_type("M", {"W"}, {}, {});

    // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const double g = 1.22074408460575947536;
    const vec3r_t inc{1.0/g, 1.0/(g*g), 1.0/(g*g*g)};

    int n=3*dim*dim*dim;
    vec3r_t x{0.5,0.5,0.5};
    std::unordered_map<vec3i_t,unsigned> counts;
    for(int i=0; i<n; i++){
        b.add_monomer("M", x*double(dim));
        counts[vec3_floor(x*dim)] += 1;
        
        for(int d=0; d<3; d++){
            x[d] += inc[d];
            if(x[d]>=1){
                x[d] -= 1;
            }
        }
    }

    std::vector<std::pair<unsigned,vec3i_t>> all;
    for(const auto &z : counts){
        all.push_back({z.second, z.first});
    }
    std::sort(all.begin(), all.end(), [](auto a, auto b){ return a.first < b.first; });
    for(int i=0; i<10; i++){
        std::cerr<<all[i].first<<" : "<<all[i].second<<"\n";
    }
    for(int i=0; i<10; i++){
        std::cerr<<all[all.size()-1-i].first<<" : "<<all[all.size()-1-i].second<<"\n";
    }

    return b.extract();
}

int main(int argc, const char *argv[])
{
    int max_parallelism=tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    if(getenv("PDPD_NUM_THREADS")){
	    max_parallelism=std::atoi(getenv("PDPD_NUM_THREADS"));
      }
      std::cerr<<"TBB is using "<<max_parallelism<<" threads.\n";
      tbb::global_control tbb_control_threads(tbb::global_control::max_allowed_parallelism, max_parallelism);
      

    std::string engine_name;
    if(argc>1){
        engine_name=argv[1];
    }else{
        usage();
    }

    std::string mode;
    if(argc>2){
        mode=argv[2];
    }else{
        usage();
    }

    int nStart=1;
    if(argc>3){
        nStart=std::stoi(argv[3]);
    }

    WorldState state;

    if(mode.substr(0,7)=="uniform"){
        auto nstr=mode.substr(8);
        int n=stoi(nstr);
        fprintf(stderr, "n=%u\n", n);

        state=make_uniform(n);
    }else{
        fprintf(stderr, "Unknown mode %s\n", mode.c_str());
    }

    std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);

    int volume=state.box[0]*state.box[1]*state.box[2];


    //tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 1);

    double tRatePrev=0;
    for(int todo=nStart; todo<1000000; todo=std::max(2,todo*3/2)){
        double t0=now();

        engine->Attach(&state);

        double t1=now();

        engine->Run(todo);

        double t2=now();

        DPDEngine::timings_t timings;
        bool timings_valid= engine->GetTimings(timings);

        fprintf(stdout, "%s,%s,%d,%d,%d,%g,%g,%g",
            engine_name.c_str(), mode.c_str(), volume, (int)state.beads.size(), todo,
            t1-t0, t2-t1, state.beads.size()/(t2-t1)*todo
        );
        if(timings_valid){
            fprintf(stdout, ", %g, %g", timings.execute_to_first_bead, state.beads.size() / timings.execute_to_first_bead * todo);
        }
        fprintf(stdout, "\n");

	fflush(stdout);        

	double tRate=(t2-t0)/todo;
	fprintf(stderr, "%g\n", abs((tRate-tRatePrev)/tRatePrev) );
	if( tRatePrev>0 && (t2-t0)>10 && abs((tRate-tRatePrev)/tRatePrev) < 0.05  ){
		break;
	}
	tRatePrev=tRate;
    }

    return 0;
}
