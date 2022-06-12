
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_to_pov.hpp"

#include <fstream>

void usage()
{
    fprintf(stderr, "world_state_to_pov : input-file output-file \n");
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
        if(argc!=3){
            usage();
        }

        std::string src1_file=argv[1];
        std::string src2_file=argv[2];

        std::cerr<<"src1_file="<<src1_file<<", src2_file="<<src2_file<<"\n";

        WorldState state1=read_world_state(src1_file);

        write_to_pov(src2_file, state1);

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
