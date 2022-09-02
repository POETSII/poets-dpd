#ifndef dpd_state_to_pov_hpp
#define dpd_state_to_pov_hpp

#include "dpd_state.hpp"

#include "dpd/core/with_optional_gzip_stream.hpp"

#include <iostream>
#include <fstream>

std::ostream &write_to_pov(std::ostream &dst, const WorldState &state1)
{
    int water_bead=-1;
    for(const BeadType &bt : state1.bead_types){
        if(bt.name=="W"){
            water_bead=bt.id;
        }
    }

    auto filter=[&](const Bead &b)
    {
        return b.bead_type!=water_bead;
    };

    std::vector<const Bead *> beads;
    for(const auto &b : state1.beads){
        if(filter(b)){
            beads.push_back(&b);
        }
    }

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

        for(const auto &b : state1.beads){
            if(b.bead_type==ignored_bead_type){
                continue;
            }

            counts[b.bead_type]++;
            dst<<"sphere { < "<<b.x[0]<<", "<<b.x[1]<<", "<<b.x[2]<<" >, 0.5 texture { pigment {color "<<colour_mapping[b.bead_type]<<"} } }\n";
        }

    return dst;
}

void write_to_pov(std::string dst_path, const WorldState &state1)
{
    with_optional_gzip_ostream(dst_path, [&](std::ostream &dst){
        write_to_pov(dst, state1);
    });
}

#endif
