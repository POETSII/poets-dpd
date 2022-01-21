
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_validator.hpp"


#include <fstream>

void usage()
{
    fprintf(stderr, "world_state_to_binary : [src-file] [dst-file] \n");
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
        bool zero_forces=false;

        std::string src_file="-";
        if(argc>1){
            src_file=argv[1];
        }

        std::string dst_file="-";
        if(argc>2){
            dst_file=argv[2];
        }

        if(argc>3){
            zero_forces=atoi(argv[3]);
        }

        std::cerr<<"src_file="<<src_file<<", dst_file="<<dst_file<<", zero_forces="<<zero_forces<<"\n";

        WorldState state=read_world_state(src_file);

        validate(state, 100);

        if(zero_forces){
            for(Bead &b : state.beads){
                b.f.clear();
            }
        }

        write_world_state(dst_file, state, true);

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
