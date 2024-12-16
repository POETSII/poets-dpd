
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_to_vtk.hpp"
#include "dpd/core/dpd_state_validator.hpp"

#include "dpd/core/logging.hpp"
#include "dpd/core/logging_impl.hpp"

#include <random>
#include <fstream>

#define TBB_PREVIEW_GLOBAL_CONTROL 1

#include <tbb/global_control.h>


void usage()
{
    fprintf(stderr, "run_world : engine-name src-file output-file steps [max_r_tol]\n");
    fprintf(stderr, "  max_r_tol : Keep running batches until all bonds meet r tolerance. By default this is 0.9\n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr,"  env:\n");
    fprintf(stderr,"     PDPD_NUM_THREADS=n : Limit TBB to using n threads.\n");
    exit(1);
}

void print_exception(const std::exception& e, int level =  0)
{
    std::cerr << std::string(level, ' ') << "exception: " << e.what() << '\n';
    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e) {
        print_exception(e, level+1);
    } catch(...) {}
}
 

double now()
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

struct world_stats_t
{
    double temperature;
    double max_r;
    double max_v;
    double max_f;
};

int main(int argc, const char *argv[])
{
    try{
        double tStart=now();

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

        std::string srcFile;
        if(argc>2){
            srcFile=argv[2];
        }else{
            usage();
        }

        std::string dstFile;
        if(argc>3){
            dstFile=argv[3];
        }else{
            usage();
        }

        int steps=10;
        if(argc>4){
            steps=std::stoi(argv[4]);
        }

        double max_r_tol=0.9;
        if(argc>5){
            max_r_tol=strtod(argv[5], 0);
            if(max_r_tol < 0.1 || max_r_tol >= 2){
                fprintf(stderr, "max_r_tol = %f doesn't make sense (?)\n", max_r_tol);
                exit(1);
            }
        }

        std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);
        engine->PrepareResources();

        double tBegin=now();
        WorldState state=read_world_state(srcFile);

        validate(state, engine->GetMaxBondLength());

        for(auto &pt : state.polymer_types){
            for(auto &b : pt.bonds){
                if(b.r0 > max_r_tol*0.9){
                    fprintf(stderr, "max_r_tol is %f, but bond pair has r0 = %f. Not going to hit stable state.", b.r0, max_r_tol);
                    exit(1);
                }
            }
        }

        std::string ok=engine->CanSupport(&state);
        if(!ok.empty()){
            fprintf(stderr, "Engine can't support world : %s\n", ok.c_str());
            exit(1);
        }

        engine->Attach(&state);

        // Max speed allowed is 0.8 box per time-step
        double max_v_limit = 0.8 / state.dt;
        // Max force would cause speed to increase by max_v in one time-step
        double max_f_limit = max_v_limit / state.dt;
        // or distance to change by more than 0.8 in one time-step
        max_f_limit = std::min(max_f_limit, 1.6/(state.dt*state.dt));

        auto calc_worst = [&]() -> world_stats_t
        {

            double sumSqMom=0;
            for(const auto &b : state.beads){
                for(int d=0; d<3; d++){
                    sumSqMom += b.x[d] * b.x[d];
                }
            }
            double temperature=sumSqMom/(3.0*state.beads.size());

            double worst_r=0;
            double sum_r=0;
            double n_r=0;
            double max_v=0;
            double max_f=0;
            for(const auto &p : state.polymers){
                const auto &pt=state.polymer_types[p.polymer_type];
                for(const auto &b : pt.bonds){
                    const auto &head=state.beads[p.bead_ids[b.bead_offset_head]];
                    const auto &tail=state.beads[p.bead_ids[b.bead_offset_tail]];

                    double r=vec3_wrapped_distance(head.x, tail.x, state.box);
                    worst_r=std::max(worst_r, r);
                    sum_r += r;
                    n_r += 1;
                }

                for(auto id : p.bead_ids){
                    const auto &b = state.beads[id];
                    max_v = std::max(max_v, b.v.l2_norm());
                    max_f = std::max(max_f, b.v.l2_norm());
                }
            }
            world_stats_t stats;
            stats.temperature=temperature;
            stats.max_r = worst_r;
            stats.max_v = max_v;
            stats.max_f = max_f;

            fprintf(stderr, "  t=%d, Temperature=%g, Angles: Target=%f, worst=%f, mean=%f, Velocity: max=%f, limit=%f, Force: max=%f, limit=%f\n", state.t, temperature, max_r_tol, worst_r, sum_r/n_r, stats.max_v, max_v_limit, stats.max_f, max_f_limit);
            return stats;
        };

        calc_worst();

        unsigned batches=0;
        while(1){
            ++batches;
            fprintf(stderr, "Starting batch %d\n", batches);
            engine->Run(steps);

            auto state=calc_worst();
            if( (state.max_r < max_r_tol)
                &&
                (state.max_v < max_v_limit)
                &&
                (state.max_f < max_f_limit)
            ){
                break;
            }
        }

        write_world_state(dstFile, state);

        fprintf(stderr, "Successful completion.\n");

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
