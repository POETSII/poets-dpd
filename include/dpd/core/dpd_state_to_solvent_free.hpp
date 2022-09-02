#ifndef dpd_core_dpd_state_to_solvent_free_hpp
#define dpd_core_dpd_state_to_solvent_free_hpp

#include <iostream>

#include "dpd/core/dpd_state.hpp"
#include "dpd/core/with_optional_gzip_stream.hpp"

std::ostream &write_to_solvent_free(std::ostream &dst, const WorldState &s)
{
    int water_bead=-1;
    for(const BeadType &bt : s.bead_types){
        if(bt.name=="W"){
            water_bead=bt.id;
        }
    }

    for(const Bead & b : s.beads){
        if(b.bead_type==water_bead){
            continue;
        }

        dst<<(int)(b.polymer_id+1)<<" "<<(int)b.polymer_type<<" "<<(int)(b.bead_id+1)<<" "<<(int)b.bead_type;
        dst<<" "<<s.bead_types[b.bead_type].r;
        dst<<" "<<b.x[0]<<" "<<b.x[1]<<" "<<b.x[2];
        dst<<" "<<b.x[0]<<" "<<b.x[1]<<" "<<b.x[2];
        dst<<" "<<b.v[0]<<" "<<b.v[1]<<" "<<b.v[2];
        dst<<"\n";
    }
    return dst;
}

void write_to_solvent_free(const std::string &dst, const WorldState &s)
{
    with_optional_gzip_ostream(dst, [&](std::ostream &dst_stream){
        write_to_solvent_free(dst_stream, s);
    });
}


#endif
