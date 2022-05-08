
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
    fprintf(stderr, "run_world : engine-name src-file output-base-name interval_count state_interval_size [vtk_interval_size] \n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr,"  env:\n");
    fprintf(stderr,"     PDPD_LOG=log-path : do full force logging to given file.\n");
    fprintf(stderr,"     PDPD_NUM_THREADS=n : Limit TBB to using n threads.\n");
    fprintf(stderr,"\n");
    fprintf(stderr, "    Either vtk_interval_size must divide state_interval_size, or vice-versa\n");
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

        std::string baseName;
        if(argc>3){
            baseName=argv[3];
        }else{
            usage();
        }

        std::string log_dst;
        if(getenv("PDPD_LOG")){
            log_dst=getenv("PDPD_LOG");
            ForceLogging::set_logger(new FileLogger(log_dst));
            fprintf(stderr, "Logging to %s\n", log_dst.c_str());
        }

        int interval_count=1000;
        if(argc>4){
            interval_count=std::stoi(argv[4]);
        }

        int state_interval_size=1;
        if(argc>5){
            state_interval_size=std::stoi(argv[5]);
        }
   
        int vtk_interval_size=state_interval_size;
        if(argc>6){
            vtk_interval_size=std::atoi(argv[6]);
        }

        int interval_size=std::min(state_interval_size,vtk_interval_size);

        std::cerr<<"interval_count="<<interval_count<<", interval_size="<<interval_size<<", state_interval_size="<<state_interval_size<<", vtk_interval_size="<<vtk_interval_size<<"\n";

        if( ! (((state_interval_size%vtk_interval_size)==0) || ((vtk_interval_size%state_interval_size)==0)) ){
           fprintf(stderr, "state_interval_size and vtk_interval_size are not compatible.\n");
           exit(1);
        }

        WorldState state=read_world_state(srcFile);

        std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);

        validate(state, engine->GetMaxBondLength());

        std::string ok=engine->CanSupport(&state);
        if(!ok.empty()){
            fprintf(stderr, "Engine can't support world : %s\n", ok.c_str());
            exit(1);
        }

        int volume=state.box[0]*state.box[1]*state.box[2];

        double t0=now();

        engine->Attach(&state);

        double t1=now();

        int state_modulus=state_interval_size / interval_size;
        int vtk_modulus=vtk_interval_size / interval_size;

        unsigned done=0;
        unsigned slice_i=0;
        engine->Run(interval_count, interval_size, [&](){
	    ++slice_i;

            double t2=now();
            done+=interval_size;

            fprintf(stdout, "%s,%s,%d,%d,%d,%g,%g,%g,%g\n",
                engine_name.c_str(), srcFile.c_str(), volume, (int)state.beads.size(), done,
                t1-t0, t2-t1, state.beads.size()/(t2-t0)*done, state.beads.size()/(t2-t1)*interval_size
            );
            fflush(stdout);
            t1=t2;

            std::vector<char> tmp(baseName.size()+100);
            std::ofstream output;

            if( (slice_i%state_modulus) == 0 ){
        	    snprintf(&tmp[0], tmp.size()-1, "%s.%09d.state", baseName.c_str(), (int)state.t);
                output.open(&tmp[0]);

                    if(!output.is_open()){
                        fprintf(stderr, "Couldnt create file %s\n", &tmp[0]);
                        exit(1);
                    }
                    write_world_state(output, state);
                    output.close();
            }

            if( (slice_i%vtk_modulus) == 0){
             snprintf(&tmp[0], tmp.size()-1, "%s.%09d.vtk", baseName.c_str(), (int)state.t);
             output.open(&tmp[0]);
             if(!output.is_open()){
                 fprintf(stderr, "Couldnt create file %s\n", &tmp[0]);
                 exit(1);
             }
             write_to_vtk(output, state);
             output.close();
	   }



            return true;
        });

        if(ForceLogging::logger()){
            delete ForceLogging::logger();
            ForceLogging::logger()=0;
        }

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
