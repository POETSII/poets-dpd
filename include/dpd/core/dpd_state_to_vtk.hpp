#ifndef dpd_state_to_vtk_hpp
#define dpd_state_to_vtk_hpp

#include "dpd_state.hpp"

#include <iostream>
#include <fstream>

std::ostream &write_to_vtk(std::ostream &dst, const WorldState &s)
{
    int water_bead=0;
    for(const BeadType &bt : s.bead_types){
        if(bt.name=="W"){
            water_bead=bt.id;
        }
    }

    auto filter=[&](const Bead &b)
    {
        return b.bead_type!=water_bead;
    };

    std::vector<const Bead *> beads;
    for(const auto &b : s.beads){
        if(filter(b)){
            beads.push_back(&b);
        }
    }

    dst<<"# vtk DataFile Version 2.0\n";
	dst<<"data\n";
	dst<<"ASCII\n";
	dst<<"DATASET POLYDATA\n";
	dst<<"POINTS "<<beads.size()<<" double\n";
    for(unsigned i=0; i<beads.size(); i++){
        const auto &b=*beads[i];
        dst<<b.x[0]<<" "<<b.x[1]<<" "<<b.x[2]<<"\n";
    }

    dst<<"POINT_DATA "<<beads.size()<<"\n";

    dst<<"SCALARS bead_type double\n";
    dst<<"LOOKUP_TABLE default\n";
    for(unsigned i=0; i<beads.size(); i++){
        dst<<(int)(beads[i]->bead_type)<<"\n";
    }

/*
    
    dst<<"VECTORS v double\n";
    for(unsigned i=0; i<beads.size(); i++){
        const auto &b=*beads[i];
        dst<<b.v[0]<<" "<<b.v[1]<<" "<<b.v[2]<<"\n";
    }

    dst<<"VECTORS f double\n";
    for(unsigned i=0; i<beads.size(); i++){
        const auto &b=*beads[i];
        dst<<b.f[0]<<" "<<b.f[1]<<" "<<b.f[2]<<"\n";
    }
*/

    return dst;
}

struct VTKSnapshotter
{
    std::string base_name;
    int sequence_number=0;
    double next_snapshot=0;
    double delta_snapshot=0;

    VTKSnapshotter(const std::string &_base_name, double _dt)
    {
        base_name=_base_name;
        delta_snapshot=_dt;
    }

    void capture(const WorldState &s)
    {
        if(s.t>=next_snapshot){
            std::string name=base_name+"."+std::to_string(sequence_number);
            std::ofstream dst(name);
            if(!dst.is_open()){
                throw std::runtime_error("Couldn't open file "+name);
            }

            write_to_vtk(dst, s);

            next_snapshot += delta_snapshot;
            sequence_number++;
        }
    }
};

#endif
