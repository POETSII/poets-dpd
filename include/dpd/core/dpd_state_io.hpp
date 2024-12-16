#ifndef dpd_state_io_hpp
#define dpd_state_io_hpp

#include "dpd_state.hpp"
#include "split_string.hpp"
#include "with_optional_gzip_stream.hpp"

#include "with_optional_gzip_stream.hpp"

#include "cvec3.hpp"

#include <iostream>

std::ostream &write_bead_type(std::ostream &dst, const BeadType &b, const WorldState &)
{
    dst<<"BeadType "<<b.id<<" "<<b.name<<" "<<b.r;
    if(b.stationary){
        dst<<" stationary";
    }
    dst<<"\n";
    return dst;
}

BeadType read_bead_type(std::istream &src, int &line_no)
{
    auto p=read_prefixed_line_and_split_on_space(src, "BeadType", -1, line_no);
    if(p.parts.size()!=4 && p.parts.size()!=5){
        throw std::runtime_error("Expecting 4 or 5 elements for prefix BeadType on line "+std::to_string(line_no));
    }
    BeadType res;
    res.id=p.unsigned_at(1);
    res.name=p.string_at(2);
    res.r=p.double_at(3);
    res.stationary=false;
    if(p.parts.size()>4){
        std::string s=p.string_at(4);
        if(s=="stationary"){
            res.stationary=true;
        }else{
            throw std::runtime_error("Couldn't intepret BeadType flag "+s+" at line "+std::to_string(line_no));
        }
    }
    return res;
}

std::ostream &write_bead(std::ostream &dst, const Bead &b, const WorldState &s)
{
    dst<<"B "<<b.bead_id<<" "<<b.polymer_id<<" "<<s.polymers.at(b.polymer_id).polymer_type<<" "<<(unsigned)b.polymer_offset;
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
        if(res.x[i]==s.box[i]){
            res.x[i]=0; // It may have been rounded up due to output rounding.
        }
    }
    for(unsigned i=0; i<3; i++){
        res.v[i]=p.double_at(i+8);
    }
    for(unsigned i=0; i<3; i++){
        res.f[i]=p.double_at(i+11);
    }

    s.beads[res.bead_id]=res;
}

///////////////////////////////////////////////////////////////////////////////////////////
// V1+V2 binary

/*
    Positions in binary are encoded as 32-bit numbers with 16-bits of fractional precision.
    They are stored relative to the origin, the same as WorldState, but not like the text format.

     No attempt is made to delta-encode or anything, though the MSBs should mostly be zero.
    
    Velocity and force are encoded at 8 bytes per vector using CVec3.

    Typically this results in 2 + 12 + 8 + 8 = 30 bytes per bead, versus ~100 bytes per bead for the text version with default cout precision

    For v2 we encode force and velocity as CVec3Half, so about 22 bytes per bead.

    Leading zeros should be common in the header and position, so can probably be compressed a little with conventional methods.
*/

struct BinaryIOWriteContext
{
    std::ostream &dst;
    uint32_t prev_polymer_id=0xFFFFFFFFul;
    uint32_t prev_bead_id=0xFFFFFFFFul;

    bool use_v2_binary=false;

    template<class T>
    void write(const T &x)
    {
        dst.write((char*)&x, sizeof(T));
        if(dst.bad()){
            throw std::runtime_error("Error while writing binary stream.");
        }
    }

    void write_u32(uint32_t x)
    {
        write(x);
    }
};

struct PolymerBinary
{
    uint16_t polymer_type : 13;
    uint16_t non_sequential_polymer_id : 1;    // If 0 then the polymer_id is previous + 1
    uint16_t non_sequential_bead_id_inter : 1; // If 0 then first bead id is previous bead id + 1
    uint16_t non_sequential_bead_id_intra : 1; // If 0 then bead ids in polymer are sequential
};

struct BinaryIOReadContext
{
    std::istream &src;
    uint32_t prev_polymer_id=0xFFFFFFFFul;
    uint32_t prev_bead_id=0xFFFFFFFFul;

    bool use_v2_binary=false;

    template<class T>
    void read(T &x)
    {
        src.read((char*)&x, sizeof(T));
        if(src.bad()){
            throw std::runtime_error("Error while reading binary stream.");
        }
    }

    uint32_t read_u32()
    {
        uint32_t x;
        read(x);
        return x;
    }
};

void write_polymer_binary(const Polymer &p, const WorldState &s, BinaryIOWriteContext &ctxt)
{
    PolymerBinary header;

    header.non_sequential_polymer_id=p.polymer_id != ctxt.prev_polymer_id+1;
    header.non_sequential_bead_id_inter=p.bead_ids.at(0) != ctxt.prev_bead_id+1;
    header.non_sequential_bead_id_intra=false;
    for(unsigned i=1; i<p.bead_ids.size(); i++){
        if( p.bead_ids[i-1]+1 != p.bead_ids[i] ){
            header.non_sequential_bead_id_intra=true;
            break;
        }
    }
    header.polymer_type=p.polymer_type;
        
    ctxt.write(header);
    if(header.non_sequential_polymer_id){
        ctxt.write_u32(p.polymer_id-ctxt.prev_polymer_id);
    }
    ctxt.prev_polymer_id=p.polymer_id;

    const double POS_SCALE=65536; // Use 16-bits of fractional precision

    bool first=true;
    for(unsigned bead_id : p.bead_ids){
        if( (first && header.non_sequential_bead_id_inter) || (!first && header.non_sequential_bead_id_intra)){
            ctxt.write_u32(bead_id - ctxt.prev_bead_id);
        }
        ctxt.prev_bead_id=bead_id;

        const auto &bead=s.beads.at(bead_id);

        std::array<uint32_t,3> pos;
        for(unsigned d=0; d<3; d++){
            double rr=round(bead.x[d]*POS_SCALE);
            if( rr == s.box[d]*POS_SCALE ){
                rr=0;
            }
            if(rr < 0 || rr == s.box[d]*POS_SCALE){
                throw std::runtime_error("Bead is out of bounds.");
            }
            pos[d]=(uint32_t)rr;
        };

        ctxt.write(pos);

        if(ctxt.use_v2_binary){
            CVec3Half v(bead.v);
            ctxt.write(v);

            CVec3Half f(bead.f);
            ctxt.write(f);

        }else{
            CVec3 v(bead.v);
            ctxt.write(v);

            CVec3 f(bead.f);
            ctxt.write(f);
        }

        first=false;
    }

    ctxt.prev_bead_id=p.bead_ids.back();
}

void read_polymer_binary(WorldState &s, BinaryIOReadContext &ctxt)
{
    PolymerBinary header;
    ctxt.read(header);

    uint32_t polymer_id;
    if(header.non_sequential_polymer_id){
        std::cerr<<"Non seq polymer id\n";
        polymer_id = ctxt.prev_polymer_id + ctxt.read_u32();
    }else{
        polymer_id = ctxt.prev_polymer_id + 1;
    }
    ctxt.prev_polymer_id = polymer_id;
    //std::cerr<<"Begin polymer "<<polymer_id<<"\n";

    if(polymer_id >= s.polymers.size()){
        throw std::runtime_error("Polymer id out of range.");
    }

    auto &poly=s.polymers.at(polymer_id);
    if(poly.bead_ids.size()!=0){
        throw std::runtime_error("Polymer appears twice.");
    }
    poly.polymer_id=polymer_id;
    poly.polymer_type=header.polymer_type;
    const auto &pt=s.polymer_types.at(poly.polymer_type);

    bool first=true;
    for(unsigned i=0; i<pt.bead_types.size(); i++){
        uint32_t bead_id;
        if( (first && header.non_sequential_bead_id_inter) || (!first && header.non_sequential_bead_id_intra) ){
            bead_id=ctxt.prev_bead_id + ctxt.read_u32();
        }else{
            bead_id=ctxt.prev_bead_id+1;
        }
        ctxt.prev_bead_id=bead_id;

        poly.bead_ids.push_back(bead_id);

        Bead &bead=s.beads.at(bead_id);

        if(bead.bead_id != (uint32_t)-1){
            throw std::runtime_error("Bead appears twice.");
        }
        bead.bead_id=bead_id;
        bead.bead_type=pt.bead_types[i];
        bead.is_monomer=pt.bead_types.size()==1;
        bead.polymer_id=polymer_id;
        bead.polymer_offset=i;
        bead.polymer_type=poly.polymer_type;

        std::array<uint32_t,3> pos;
        ctxt.read(pos);
        for(int d=0; d<3; d++){
            double x=ldexp(pos[d],-16);
            if(x < 0 || x >= s.box[d]){
                throw std::runtime_error("Bead was out of box.");
            }

            bead.x[d]=x;
        }

        if(ctxt.use_v2_binary){
            CVec3Half v;
            ctxt.read(v);
            bead.v=v.get_vec3r();

            CVec3Half f;
            ctxt.read(f);
            bead.f=f.get_vec3r();
        }else{
            CVec3 v;
            ctxt.read(v);
            bead.v=v.get_vec3r();

            CVec3 f;
            ctxt.read(f);
            bead.f=f.get_vec3r();
        }

        first=false;
    }
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
            res.bonds[i].kappa=p.double_at(3);
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
            res.bond_pairs[i].kappa=p.double_at(3);
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

std::ostream &write_world_state(std::ostream &dst, const WorldState &state, bool binary=false)
{
    bool use_v2_binary=true;
        
    if(binary){
        dst<<"WorldState v2binary "<<state.bead_types.size()<<" "<<state.polymer_types.size()<<" "<<state.beads.size()<<" "<<state.polymers.size()<<"\n";
        //dst<<"WorldState v1binary "<<state.bead_types.size()<<" "<<state.polymer_types.size()<<" "<<state.beads.size()<<" "<<state.polymers.size()<<"\n";
    }else{
        dst<<"WorldState v0 "<<state.bead_types.size()<<" "<<state.polymer_types.size()<<" "<<state.beads.size()<<" "<<state.polymers.size()<<"\n";
    }
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

    if(binary){
        dst<<"BinaryData\n";
        char ch=0;
        dst.write(&ch, 1);
        BinaryIOWriteContext ctxt{dst};
        ctxt.use_v2_binary=use_v2_binary;
        for(const Polymer &p : state.polymers){
            write_polymer_binary(p, state, ctxt);
        }
    }else{
        dst<<"# Beads\n";
        for(const auto &b : state.beads){
            write_bead(dst, b, state);
        }
        dst<<"\n";
    }
    return dst;
}

void write_world_state(std::string dst, const WorldState &state, bool binary=false)
{
    with_optional_gzip_ostream(dst, [&](std::ostream &dst){
        write_world_state(dst, state, binary);
    });
}

WorldState read_world_state(std::istream &src, int &line_no)
{
    WorldState res;

    auto p=read_prefixed_line_and_split_on_space(src, "WorldState", 6, line_no);
    int world_state_line_no=line_no;
    try{
        bool is_binary=false;
        bool use_v2_binary=false;
        if(p.string_at(1)=="v2binary"){
            is_binary=true;
            use_v2_binary=true;
        }else if(p.string_at(1)=="v1binary"){
            is_binary=true;
        }else if(p.string_at(1)!="v0"){
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
                res.interactions[i*numBeadTypes+j].conservative=p.expression_at(2+j);
            }
        }
        for(unsigned i=0; i<numBeadTypes; i++){
            p=read_prefixed_line_and_split_on_space(src, "DissipativeStrength", numBeadTypes+2, line_no);
            for(unsigned j=0; j<numBeadTypes; j++){
                res.interactions[i*numBeadTypes+j].dissipative=p.expression_at(2+j);
            }
        }

        for(auto &bt : res.bead_types){
            bt=read_bead_type(src, line_no);
        }
        for(auto & s: res.polymer_types){
            read_polymer_type(src, line_no, res);
        }

        if(is_binary){
            auto p=read_prefixed_line_and_split_on_space(src, "BinaryData", 1, line_no);
            while(1){
                char ch;
                src.read(&ch, 1);
                if(src.bad()){
                    throw std::runtime_error("Missing null before binary data.");
                }
                if(ch==0){
                    break;
                }
                if(!isspace(ch)){
                    throw std::runtime_error("Unexpected char before null char.");
                }                
            }

            BinaryIOReadContext ctxt{src};
            ctxt.use_v2_binary=use_v2_binary;
            for(unsigned i=0; i<res.polymers.size(); i++){
                read_polymer_binary(res, ctxt);
            }
        }else{
            for(unsigned i=0; i<res.beads.size(); i++){
                read_bead(src, line_no, res);
            }
        }

    }catch(...){
        std::throw_with_nested(std::runtime_error("Exception while parsing WorldState at line '"+std::to_string(world_state_line_no)));
    }

    return res;
}

WorldState read_world_state(std::string src)
{
    WorldState res;
    with_optional_gzip_istream(src, [&](std::istream &src){
        int line_no=0;
        res=read_world_state(src, line_no);
    });
    return res;
}

#endif
