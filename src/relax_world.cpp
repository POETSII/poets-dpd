
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

        auto calc_worst = [&]() -> double
        {

            double sumSqMom=0;
            for(const auto &b : state.beads){
                for(int d=0; d<3; d++){
                    sumSqMom += b.x[d] * b.x[d];
                }
            }
            double temperature=sumSqMom/(3.0*state.beads.size());

            double worst=0;
            double sum=0;
            double n=0;
            for(const auto &p : state.polymers){
                const auto &pt=state.polymer_types[p.polymer_type];
                for(const auto &b : pt.bonds){
                    const auto &head=state.beads[p.bead_ids[b.bead_offset_head]];
                    const auto &tail=state.beads[p.bead_ids[b.bead_offset_tail]];

                    double r=vec3_wrapped_distance(head.x, tail.x, state.box);
                    worst=std::max(worst, r);
                    sum += r;
                    n += 1;
                }
            }
            fprintf(stderr, "  t=%d, Temperature=%g, Angles: Target=%f, worst=%f, mean=%f\n", state.t, temperature, max_r_tol, worst, sum/n);
            return worst;
        };

        calc_worst();

        unsigned batches=0;
        while(1){
            ++batches;
            fprintf(stderr, "Starting batch %d\n", batches);
            engine->Run(steps);

            double worst=calc_worst();
            if(worst < max_r_tol){
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
