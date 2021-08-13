#ifndef dpd_state_io_hpp
#define dpd_state_io_hpp

#include "dpd_state.hpp"
#include "split_string.hpp"

#include <iostream>

std::ostream &write_bead_type(std::ostream &dst, const BeadType &b, const WorldState &)
{
    return dst<<"BeadType "<<b.id<<" "<<b.name<<" "<<b.r<<"\n";
}

BeadType read_bead_type(std::istream &src, int &line_no)
{
    auto p=read_prefixed_line_and_split_on_space(src, "BeadType", 4, line_no);
    BeadType res;
    res.id=p.unsigned_at(1);
    res.name=p.string_at(2);
    res.r=p.double_at(3);
    return res;
}

std::ostream &write_bead(std::ostream &dst, const Bead &b, const WorldState &s)
{
    dst<<"B "<<b.bead_id<<" "<<b.polymer_id<<" "<<s.polymers.at(b.polymer_id).polymer_type<<" "<<b.polymer_offset;
    dst<<" ";
    for(unsigned i=0; i<3; i++){
        dst<<" "<<b.x[i]+s.origin[i];
    }
    dst<<" ";
    for(unsigned i=0; i<3; i++){
        dst<<" "<<b.v[i];
    }
    dst<<" ";  
    for(unsigned i=0; i<3; i++){
        dst<<" "<<b.f[i];
    }
    dst<<"\n";
    return dst;
}

void read_bead(std::istream &src, int line_no, WorldState &s)
{
    Bead res;
    auto p=read_prefixed_line_and_split_on_space(src, "B", 14, line_no);
    
    res.bead_id=p.unsigned_at(1);
    res.polymer_id=p.unsigned_at(2);
    res.polymer_type=p.unsigned_at(3);
    res.polymer_offset=p.unsigned_at(4);

    const auto &poly_type=s.polymer_types.at(res.polymer_type); 
    auto &poly=s.polymers.at(res.polymer_id);
    if(poly.polymer_id==(uint32_t)-1){
        poly.polymer_id=res.polymer_id;
        poly.polymer_type=res.polymer_type;
        poly.bead_ids.assign(poly_type.bead_types.size(), -1);
    }

    res.is_monomer=poly_type.bead_types.size()==1;

    res.bead_type=poly_type.bead_types.at(res.polymer_offset);
    
    if(poly.bead_ids.at(res.polymer_offset) != (uint32_t)-1 ){
        throw std::runtime_error("Duplicate bead within polymer");
    }
    poly.bead_ids.at(res.polymer_offset) = res.bead_id;

    for(unsigned i=0; i<3; i++){
        res.x[i]=p.double_at(i+5) - s.origin[i];
    }
    for(unsigned i=0; i<3; i++){
        res.v[i]=p.double_at(i+8);
    }
    for(unsigned i=0; i<3; i++){
        res.f[i]=p.double_at(i+11);
    }

    s.beads[res.bead_id]=res;
}

std::ostream &write_polymer_type(std::ostream &dst, const PolymerType &m, const WorldState &)
{
    dst<<"PolymerType "<<m.polymer_id<<" "<<m.name<<" "<<m.bead_types.size()<<" "<<m.bonds.size()<<" "<<m.bond_pairs.size()<<"\n";
    dst<<"BeadTypeIndices";
    for(const auto &bt : m.bead_types){
        dst<<" "<<bt;
    }
    dst<<"\n";
    for(const auto &b : m.bonds){
        dst<<"Bond "<<b.bead_offset_head<<" "<<b.bead_offset_tail<<" "<<b.kappa<<" "<<b.r0<<"\n";
    }
    for(const auto &bp : m.bond_pairs){
        dst<<"BondPair "<<bp.bond_offset_head<<" "<<bp.bond_offset_tail<<" "<<bp.kappa<<" "<<bp.theta0<<"\n";
    }
    dst<<"\n";
    return dst;
}


void read_polymer_type(std::istream &src, int &line_no, WorldState &state)
{
    auto p=read_prefixed_line_and_split_on_space(src, "PolymerType", 6, line_no);
    int polymer_line_no=line_no;

    try{
        int polymer_id=p.unsigned_at(1);

        PolymerType &res=state.polymer_types.at(polymer_id);
        if(res.polymer_id!=(uint32_t)-1){
            throw std::runtime_error("Duplicate polymer.");
        }

        res.polymer_id=polymer_id;
        res.name=p.string_at(2);
        res.bead_types.resize(p.unsigned_at(3));
        res.bonds.resize(p.unsigned_at(4));
        res.bond_pairs.resize(p.unsigned_at(5));

        p=read_prefixed_line_and_split_on_space(src, "BeadTypeIndices", 1+res.bead_types.size(), line_no);
        for(unsigned i=0; i<res.bead_types.size(); i++){
            res.bead_types[i]=p.unsigned_at(i+1);
        }
        for(unsigned i=0; i<res.bonds.size(); i++){
            p=read_prefixed_line_and_split_on_space(src, "Bond", 5, line_no);
            //dst<<"Bond "<<b.bead_offset_a<<" "<<b.bead_offset_b<<" "<<b.kappa<<" "<<b.r0<<"\n";
            res.bonds[i].bead_offset_head=p.unsigned_at(1);
            res.bonds[i].bead_offset_tail=p.unsigned_at(2);
            res.bonds[i].kappa=p.unsigned_at(3);
            res.bonds[i].r0=p.double_at(4);

            if(res.bonds[i].bead_offset_head >= res.bead_types.size()){
                throw std::runtime_error("Invalid bead offset "+std::to_string(res.bonds[i].bead_offset_head) + " at line "+std::to_string(line_no));
            }
            if(res.bonds[i].bead_offset_tail >= res.bead_types.size()){
                throw std::runtime_error("Invalid bead offset b "+std::to_string(res.bonds[i].bead_offset_tail) + " at line "+std::to_string(line_no));
            }
        }
        for(unsigned i=0; i<res.bond_pairs.size(); i++){
            p=read_prefixed_line_and_split_on_space(src, "BondPair", 5, line_no);
            // dst<<"BondPair "<<bp.bond_offset_a<<" "<<bp.bond_offset_b<<" "<<bp.kappa<<" "<<bp.theta0<<"\n";
            res.bond_pairs[i].bond_offset_head=p.unsigned_at(1);
            res.bond_pairs[i].bond_offset_tail=p.unsigned_at(2);
            res.bond_pairs[i].kappa=p.unsigned_at(3);
            res.bond_pairs[i].theta0=p.double_at(4);

            if(res.bond_pairs[i].bond_offset_head >= res.bonds.size()){
                throw std::runtime_error("Invalid bond offset a "+std::to_string(res.bond_pairs[i].bond_offset_head) + " at line "+std::to_string(line_no));
            }
            if(res.bond_pairs[i].bond_offset_tail >= res.bonds.size()){
                throw std::runtime_error("Invalid bond offset b "+std::to_string(res.bond_pairs[i].bond_offset_tail) + " at line "+std::to_string(line_no));
            }
        }

    }catch(...){
        std::throw_with_nested(std::runtime_error("Exception while parsing PolymerType at line '"+std::to_string(polymer_line_no)));
    }
}

std::ostream &write_world_state(std::ostream &dst, const WorldState &state)
{
    dst<<"WorldState v0 "<<state.bead_types.size()<<" "<<state.polymer_types.size()<<" "<<state.beads.size()<<" "<<state.polymers.size()<<"\n";
    dst<<"T "<<state.t<<" "<<state.dt<<"\n";
    dst<<"Lambda "<<state.lambda<<"\n";
    dst<<"Origin "<<state.origin[0]<<" "<<state.origin[1]<<" "<<state.origin[2]<<"\n";
    dst<<"Box "<<state.box[0]<<" "<<state.box[1]<<" "<<state.box[2]<<"\n";
    dst<<"Seed "<<state.seed<<"\n";
    unsigned numBeadTypes=state.bead_types.size();
    for(unsigned i=0; i<numBeadTypes; i++){
        dst<<"ConservativeStrength "<<i;
        for(unsigned j=0; j<numBeadTypes; j++){
            dst<<" "<<state.interactions[i*numBeadTypes+j].conservative;
        }
        dst<<"\n";
    }
    for(unsigned i=0; i<numBeadTypes; i++){
        dst<<"DissipativeStrength "<<i;
        for(unsigned j=0; j<numBeadTypes; j++){
            dst<<" "<<state.interactions[i*numBeadTypes+j].dissipative;
        }
        dst<<"\n";
    }
    dst<<"\n";

    dst<<"# BeadTypes\n";
    for(const auto &bt : state.bead_types){
        write_bead_type(dst, bt, state);
    }
    dst<<"\n";
    
    dst<<"# polymerTypes\n";
    for(const auto &mt : state.polymer_types){
        write_polymer_type(dst, mt, state);
    }
    dst<<"\n";

    dst<<"# Beads\n";
    for(const auto &b : state.beads){
        write_bead(dst, b, state);
    }
    dst<<"\n";
    return dst;
}

WorldState read_world_state(std::istream &src, int &line_no)
{
    WorldState res;

    auto p=read_prefixed_line_and_split_on_space(src, "WorldState", 6, line_no);
    int world_state_line_no=line_no;
    try{
        if(p.string_at(1)!="v0"){
            throw std::runtime_error("Expecting version 'v0' for WorldState, but got '"+p.string_at(1)+"'");
        }

        res.bead_types.assign(p.unsigned_at(2), {});
        res.polymer_types.assign(p.unsigned_at(3), {});
        res.beads.assign(p.unsigned_at(4), {});
        res.polymers.assign(p.unsigned_at(5), {});

        p=read_prefixed_line_and_split_on_space(src, "T", 3, line_no);
        res.t=p.double_at(1);
        res.dt=p.double_at(2);

        p=read_prefixed_line_and_split_on_space(src, "Lambda", 2, line_no);
        res.lambda=p.double_at(1);

        p=read_prefixed_line_and_split_on_space(src, "Origin", 4, line_no);
        for(unsigned i=0; i<3; i++){
            res.origin[i]=p.double_at(i+1);
        }

        p=read_prefixed_line_and_split_on_space(src, "Box", 4, line_no);
        for(unsigned i=0; i<3; i++){
            res.box[i]=p.double_at(i+1);
        }

        p=read_prefixed_line_and_split_on_space(src, "Seed", 2, line_no);
        res.seed=p.unsigned_at(1);

        unsigned numBeadTypes=res.bead_types.size();
        res.interactions.resize(numBeadTypes*numBeadTypes);
        for(unsigned i=0; i<numBeadTypes; i++){
            p=read_prefixed_line_and_split_on_space(src, "ConservativeStrength", numBeadTypes+2, line_no);
            for(unsigned j=0; j<numBeadTypes; j++){
                res.interactions[i*numBeadTypes+j].conservative=p.double_at(2+j);
            }
        }
        for(unsigned i=0; i<numBeadTypes; i++){
            p=read_prefixed_line_and_split_on_space(src, "DissipativeStrength", numBeadTypes+2, line_no);
            for(unsigned j=0; j<numBeadTypes; j++){
                res.interactions[i*numBeadTypes+j].dissipative=p.double_at(2+j);
            }
        }

        for(auto &bt : res.bead_types){
            bt=read_bead_type(src, line_no);
        }
        for(auto & s: res.polymer_types){
            read_polymer_type(src, line_no, res);
        }
        for(unsigned i=0; i<res.beads.size(); i++){
            read_bead(src, line_no, res);
        }

    }catch(...){
        std::throw_with_nested(std::runtime_error("Exception while parsing WorldState at line '"+std::to_string(world_state_line_no)));
    }

    return res;
}

#endif
