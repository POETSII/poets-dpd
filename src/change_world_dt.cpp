
#include "dpd/core/dpd_engine.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_validator.hpp"

#include <fstream>

void usage()
{
    fprintf(stderr, "change_world_dt : src-file output-file dt \n");
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

        std::string srcFile;
        if(argc>1){
            srcFile=argv[1];
        }else{
            usage();
        }

        std::string dstFile;
        if(argc>2){
            dstFile=argv[2];
        }else{
            usage();
        }

        double dt;
        if(argc>3){
            dt=std::strtod(argv[3], 0);
        }else{
		usage();
	}

        WorldState state=read_world_state(srcFile);

        validate(state, 100);

	state.dt=dt;

        write_world_state(dstFile, state);

        fprintf(stderr, "Successfully changed dt.\n");

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
