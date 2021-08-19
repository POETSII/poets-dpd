
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include <random>

void usage()
{
    fprintf(stderr, "benchmark_engine_intervals : engine-name mode interval_count interval_size \n");
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
    for(int i=0; i<n; i++){
        b.add_monomer("M", x*double(dim));
        
        for(int d=0; d<3; d++){
            x[d] += inc[d];
            if(x[d]>=1){
                x[d] -= 1;
            }
        }
    }

    return b.extract();
}

int main(int argc, const char *argv[])
{
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

    int interval_count=1;
    if(argc>3){
        interval_count=std::stoi(argv[3]);
    }

    int interval_size=1;
    if(argc>4){
        interval_size=std::stoi(argv[4]);
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

    double t0=now();

    engine->Attach(&state);

    double t1=now();

    unsigned done=0;
    engine->Run(interval_count, interval_size, [&](){
        double t2=now();
        done+=interval_size;

        fprintf(stdout, "%s,%s,%d,%d,%d,%g,%g,%g,%g\n",
            engine_name.c_str(), mode.c_str(), volume, state.beads.size(), done,
            t1-t0, t2-t1, state.beads.size()/(t2-t0)*done, state.beads.size()/(t2-t1)*interval_size
        );
        t1=t2;
        return true;
    });

    

    return 0;
}
