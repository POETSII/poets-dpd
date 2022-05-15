
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_to_vtk.hpp"

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

        std::vector<std::string> colours={ // Add more colours here if needed
            "Red", "Green", "Blue", "Yellow", "Cyan", "Magenta", "Orange", "Violet", "Wheat" 
        };
        std::vector<std::string> colour_mapping;

        int ignored_bead_type=-1;
        for(auto bt : state1.bead_types){
            if(bt.name=="W"){
                ignored_bead_type=bt.id;
                colour_mapping.push_back("Black"); // Shouldn't be used
                std::cerr<<"Ignoring bead type "<<bt.id<<", with name "<<bt.name<<"\n";
            }else{
                if(colours.empty()){
                    throw std::runtime_error("Too many bead types, not enough colours.");
                }
                colour_mapping.push_back(colours.front());
                colours.erase(colours.begin());
            }
        }

        std::cerr<<"Writing. nBeads="<<state1.beads.size()<<"\n";

        std::vector<unsigned> counts(state1.bead_types.size(), 0);

        std::ofstream dst(src2_file);

        double cx=state1.box.x[0]/2;
        double cy=state1.box.x[1]/2;
        double cz=state1.box.x[2]/2;

        double ox=state1.box.x[0]*2;
        double oy=0;
        double oz=cz;

dst<<R"(#include "colors.inc"					
#include "stones.inc"					
background {color White}
)";

        dst<<"camera{ location < "<<cx<<","<<-state1.box[1]<<","<<-cz<<" >  look_at  < "<<cx<<","<<cy<<","<<cz<<"> } \n";

        dst<<"light_source{ < "<<cx<<", "<<cy<<", "<<cz<<"> color White shadowless }\n";
        dst<<"light_source{ < "<<-cx<<", "<<0<<", "<<0<<"> color White shadowless }\n";
        dst<<"light_source{ < "<<0<<", "<<-cy<<", "<<0<<"> color White shadowless }\n";
        dst<<"light_source{ < "<<0<<", "<<0<<", "<<-cz<<"> color White shadowless }\n";
        
        dst<<"box{ < 0,0,0 > < "<<state1.box[0]<<","<<state1.box[1]<<","<<state1.box[2]<<" > texture{ pigment{ color rgbf < 0.9,0.9,0.9,0.9 > } } }\n";


//finish="""finish { ambient .2 diffuse .4 roughness .001 }"""
    std::string finish="";

        for(const auto &b : state1.beads){
            if(b.bead_type==ignored_bead_type){
                continue;
            }

            counts[b.bead_type]++;
            dst<<"sphere { <"<<b.x[0]<<", "<<b.x[1]<<", "<<b.x[2]<<" >, 0.5 texture { pigment {color "<<colour_mapping[b.bead_type]<<" } "<<finish<<" } }\n";
        }

        for(unsigned i=0; i<counts.size(); i++){
            if(i!=0){
                std::cerr<<", ";
            }
            std::cerr<<state1.bead_types[i].name<<":"<<counts[i];
        }
        std::cerr<<"\n";
    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
