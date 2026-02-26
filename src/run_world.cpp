
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_to_vtk.hpp"
#include "dpd/core/dpd_state_to_pov.hpp"
#include "dpd/core/dpd_state_to_solvent_free.hpp"
#include "dpd/core/dpd_state_validator.hpp"
#include "dpd/core/with_optional_gzip_stream.hpp"

#include "dpd/core/logging.hpp"
#include "dpd/core/logging_impl.hpp"

#include <tbb/task_group.h>

#include <random>
#include <queue>
#include <fstream>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>


void usage()
{
    fprintf(stderr, "run_world : engine-name src-file output-base-name interval_count state_interval_size [snapshot_interval_size]\n");
    fprintf(stderr, "   --vtk-snapshot : Dump all non water beads into a vtk snapshot of positions.\n");
    fprintf(stderr, "   --povray-snapshot : Dump all non water beads into a povray snapshot of positions.\n");
    fprintf(stderr, "   --solvent-free-snapshot : Dump all non water beads to an Osprey solvent-free snapshot.\n");
    fprintf(stderr, "   --povray-render : Dump positions and render to a png as well (implies --povray-snapshot).\n");
    fprintf(stderr, "   --gzip-snapshot : Snapshots (e.g. povray, vtk, solvent free) will be gzipped.\n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr,"  env:\n");
    fprintf(stderr,"     PDPD_LOG=log-path : do full force logging to given file.\n");
    fprintf(stderr,"     PDPD_NUM_THREADS=n : Limit TBB to using n threads.\n");
    fprintf(stderr,"\n");
    fprintf(stderr, "    Either snapshot_interval_size must divide state_interval_size, or vice-versa\n");
    
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

      bool vtk_snapshot=false;
      bool povray_snapshot=false;
      bool povray_render=false;
      bool solvent_free_snapshot=false;
      bool gzip_snapshot=false;
      
        std::vector<std::string> args;
        for(int i=1; i<argc; i++){
            if(!strcmp("--povray-snapshot", argv[i])){
                povray_snapshot=true;
            }else if(!strcmp("--povray-render", argv[i])){
                povray_render=true;
            }else if(!strcmp("--vtk-snapshot", argv[i])){
                vtk_snapshot=true;
            }else if(!strcmp("--solvent-free-snapshot", argv[i])){
                solvent_free_snapshot=true;
            }else if(!strcmp("--gzip-snapshot", argv[i])){
                gzip_snapshot=true;
            }else{
                args.push_back(argv[i]);
            }            
        }


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
   
        int snapshot_interval_size=state_interval_size;
        if(argc>6){
            snapshot_interval_size=std::atoi(argv[6]);
        }

        int interval_size=std::min(state_interval_size,snapshot_interval_size);

        std::cerr<<"interval_count="<<interval_count<<", interval_size="<<interval_size<<", state_interval_size="<<state_interval_size<<", snapshot_interval_size="<<snapshot_interval_size<<"\n";

        if( ! (((state_interval_size%snapshot_interval_size)==0) || ((snapshot_interval_size%state_interval_size)==0)) ){
           fprintf(stderr, "state_interval_size and snapshot_interval_size are not compatible.\n");
           exit(1);
        }

        WorldState state=read_world_state(srcFile);

        std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);

        fprintf(stderr, "Validating with max_r=%f....\n", engine->GetMaxBondLength());
        validate(state, engine->GetMaxBondLength());
        fprintf(stderr, "  ok\n");

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
        int snapshot_modulus=snapshot_interval_size / interval_size;

        std::queue<std::unique_ptr<tbb::task_group>> async_tasks;

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

            if(async_tasks.size()>4){
                async_tasks.front()->wait();
                async_tasks.pop();
            }
            async_tasks.push(std::make_unique<tbb::task_group>());

            // Create an independent shared copy for async tasks
            std::shared_ptr<WorldState> snapshot=std::make_shared<WorldState>(state);

            std::string baseNameNow;
            {
                std::vector<char> tmp(baseName.size()+100);
            	snprintf(&tmp[0], tmp.size()-1, "%s.%09d", baseName.c_str(), (int)snapshot->t);
                baseNameNow=&tmp[0];
            }

            if( (slice_i%state_modulus) == 0 ){
                async_tasks.back()->run([=](){
                    write_world_state(baseNameNow+".state.gz", *snapshot, true);
                });
            }

            if( (slice_i%snapshot_modulus) == 0){
                if(vtk_snapshot){
                    async_tasks.back()->run([=](){
                        std::string name=baseNameNow+".vtk";
                        if(gzip_snapshot){
                            name+=".gz";
                        }
                        write_to_vtk(name, *snapshot);
                    });
                }
                if(solvent_free_snapshot){
                    async_tasks.back()->run([=](){
                        std::string name=baseNameNow+".rst";
                        if(gzip_snapshot){
                            name+=".gz";
                        }
                        write_to_solvent_free(name, *snapshot);
                    });
                }
                if(povray_render || povray_snapshot){
                    async_tasks.back()->run([=,&povray_render](){
                        std::string name=baseNameNow+".pov";
                        if(gzip_snapshot){
                            name+=".gz";
                        }
                        write_to_pov(name, *snapshot);
                        if(povray_render){
                            std::string pov_name=name;
                            if(gzip_snapshot){
                                if(system( ("gunzip -k -f "+name).c_str() ) ){
                                    fprintf(stderr, "gunzip failed. Cancelling povray_render\n");
                                    povray_render=false;
                                    return;
                                }
                                pov_name=baseNameNow+".pov";
                            }
                            // Note: this ties up a thread, but is not too
                            // base for typicaly milti-core systems
                            std::stringstream cmd;
                            cmd<<"povray -D -V +O"<<baseNameNow<<".png "<<pov_name<<" 2> /dev/null";
                            if(system( cmd.str().c_str() )){
                                fprintf(stderr, "povray failed. Cancelling povray_render\n");
                                povray_render=false;
                                return;
                            }
                            if(gzip_snapshot){
                                unlink(pov_name.c_str());
                            }
                            if(!povray_snapshot){
                                unlink(name.c_str());
                            }
                        }
                    });
                }
	        }

            while(!async_tasks.empty()){
                async_tasks.front()->wait();
                async_tasks.pop();
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
