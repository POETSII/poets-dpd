
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_validator.hpp"


#include <fstream>

void usage()
{
    fprintf(stderr, "world_state_diff : src-file1 src-file1 \n");
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

        std::ifstream src1;
        src1.open(src1_file.c_str());
        if(!src1.is_open()){
            std::cerr<<"Couldn't open "<<src1_file<<"\n";
            exit(1);
        }

        std::ifstream src2;
        src2.open(src2_file.c_str());
        if(!src2.is_open()){
            std::cerr<<"Couldn't open "<<src2_file<<"\n";
            exit(1);
        }
    
        int line_no=0;
        WorldState state1=read_world_state(src1, line_no);
        validate(state1, 100);

        line_no=0;
        WorldState state2=read_world_state(src2, line_no);
        validate(state2, 100);

        #define COMPARE_PROP(p) if( ! (state1. p == state2. p )){ std::cerr<<"worlds differ on property "<< #p <<"\n"; exit(1); }

        COMPARE_PROP(box);
        COMPARE_PROP(origin);
        COMPARE_PROP(dt);
        COMPARE_PROP(lambda);
        COMPARE_PROP(interactions);
        COMPARE_PROP(lambda);
        COMPARE_PROP(seed);
        COMPARE_PROP(t);
        COMPARE_PROP(polymer_types);
        COMPARE_PROP(polymers);

        // This should be impossible if they validate and have the same polymers
        assert(state1.beads.size()==state2.beads.size());
        
        for(unsigned i=0; i<state1.beads.size(); i++){
            const auto &bead1=state1.beads[i];
            const auto &bead2=state2.beads[i];

            assert(bead1.bead_id==i);
            assert(bead2.bead_id==i);
            assert(bead1.bead_id==bead2.bead_id);

            double err= vec3_wrapped_distance(bead1.x, bead2.x, state1.box);
            if(err > 0.001){
                std::cerr<<"  Bead id "<<i<<", x1="<<bead1.x<<", x2="<<bead2.x<<", dist="<<err<<"\n";
                exit(2);
            }

            vec3r_t dist=bead1.v-bead2.v;
            err=dist.l2_norm();
            if(err > 0.001){
                std::cerr<<"  Bead id "<<i<<", v1="<<bead1.v<<", v2="<<bead2.v<<", err="<<err<<"\n";
                exit(2);
            }

            dist=bead1.f-bead2.f;
            err=dist.l2_norm();
            if(err > 0.001){
                std::cerr<<"  Bead id "<<i<<", f1="<<bead1.f<<", f2="<<bead2.f<<", err="<<err<<"\n";
                exit(2);
            }
        }


    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
