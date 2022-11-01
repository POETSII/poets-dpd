#ifndef dpd_state_to_snapshot_hpp
#define dpd_state_to_snapshot_hpp

#include "dpd_state.hpp"

#include "dpd/core/with_optional_gzip_stream.hpp"

#include <iostream>
#include <fstream>

std::ostream &write_to_snapshot(std::ostream &dst, const WorldState &state1)
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

    dst<<"{\"type\":\"dpd-snapshot-v0\", \"numSpecies\":"<<state1.bead_types.size()<<", \"numBeads\":"<<beads.size()<<",\n";
    dst<<" \"lower_bounds\":["<<state1.origin[0]<<","<<state1.origin[1]<<","<<state1.origin[2]<<"],\n";
    dst<<" \"upper_bounds\":["<<state1.origin[0]+state1.box[0]<<","<<state1.origin[1]+state1.box[1]<<","<<state1.origin[2]+state1.box[2]<<"],\n";
    dst<<" \"species\":[";
    for(int i=0; i<beads.size(); i++){
        if(i!=0)
            dst<<",";
        dst<<(int)beads[i]->bead_type;
    }
    dst<<"],\"positions\":\n[";
    for(int i=0; i<beads.size(); i++){
        if(i!=0)
            dst<<",";
        dst<<beads[i]->x[0]<<","<<beads[i]->x[1]<<","<<beads[i]->x[2];
    }
    dst<<"]\n";
    dst<<"}\n";      
}

void write_to_snapshot(const std::string &dst, const WorldState &s)
{
    with_optional_gzip_ostream(dst, [&](std::ostream &dst_stream){
        write_to_snapshot(dst_stream, s);
    });
}

#endif
