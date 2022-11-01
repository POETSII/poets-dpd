
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"

#include <fstream>

void usage()
{
    fprintf(stderr, "world_state_to_snapshot : input-file\n");
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

        std::cerr<<"src1_file="<<src1_file<<"\n";

        WorldState state1=read_world_state(src1_file);

         

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
