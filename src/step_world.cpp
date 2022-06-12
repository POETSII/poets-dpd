
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
    fprintf(stderr, "run_world : engine-name src-file [output-file|""null""] steps \n");
    fprintf(stderr, "  Writing to ""null"" will stop any export of the output state (same as /dev/null, but much faster).\n");
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

        int steps=1000;
        if(argc>4){
            steps=std::stoi(argv[4]);
        }

        std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);
        engine->PrepareResources();

        double tBegin=now();
        WorldState state=read_world_state(srcFile);

	fprintf(stderr, "Validating with max_r=%f....\n", engine->GetMaxBondLength());
        validate(state, engine->GetMaxBondLength());

        std::string ok=engine->CanSupport(&state);
        if(!ok.empty()){
            fprintf(stderr, "Engine can't support world : %s\n", ok.c_str());
            exit(1);
        }

        double load_time = now() - tBegin;

        double tStartCore=now();
        engine->Attach(&state);

	fprintf(stderr, "Stepping for %d steps at %f\n", steps, state.dt);

        engine->Run(steps);
        double tEndCore=now();
        
        if(dstFile!="null"){
            write_world_state(dstFile, state);
        }

        double tEnd=now();

        DPDEngine::timings_t timings;
        if(engine->GetTimings(timings)){
            int volume=state.box.x[0] * state.box.x[1] * state.box.x[2];
            fprintf(stdout,"# File,Engine,Steps,Volume,Beads, tLoad,tCompile,tAquire,tConfigure,tExecToFirst,tExecToLast,tPerfCounters,  bpsEndToEnd,bpsCompileConfigExecExport,bpsExec\n");
            fprintf(stdout, "%s,%s,%u,%u,%u,  ", srcFile.c_str(), engine_name.c_str(), steps, volume, state.beads.size());
            fprintf(stdout, "%f,%f,%f,%f,%f,%f,%f,  ", 
                load_time, timings.compile, timings.aquire, timings.configure, timings.execute_to_first_bead, timings.execute_to_last_bead, timings.perf_counters);
            double bpsEndToEnd= state.beads.size() * double(steps) / (tEnd-tStart);
            double bpsCompileConfigExecExport= state.beads.size() * double(steps) / (tEndCore-tStartCore);
            double bpsExec=state.beads.size() * double(steps) / timings.execute_to_first_bead;
            fprintf(stdout, "%g,%g,%g\n", bpsEndToEnd, bpsCompileConfigExecExport, bpsExec);
        }

        fprintf(stderr, "Successful completion.\n");

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
