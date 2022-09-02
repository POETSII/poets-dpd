#ifndef dpd_state_to_osprey_restart_hpp
#define dpd_state_to_osprey_restart_hpp

#include "dpd_state.hpp"
#include "dpd/core/with_optional_gzip_stream.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <mutex>

/*
    An osprey inclusive restart state is not really possible
    to recover from WorldState, as the way that polymer types
    are represented is different. Osprey wants bonds and bond-pairs
    listed in the same order as in the original dmpci, but that
    is not available in WorldState as it represents each polymer
    type as a graph.

    To get the inclusive footer we need to be given a string that
    contains the inclusive footer, which is just the bottom
    of the .dat restart file. If we have that string, the footer
    is included. Otherwise the restart state just contains beads.
*/  
std::ostream &write_to_osprey_restart(
    std::ostream &dst,
    const WorldState &s,
    const std::vector<std::string> &footer_strings={}
)
{
    // m_outStream << m_riState.IsWallPresent()    << " ";
	// m_outStream << m_riState.IsGravityPresent() << " ";
	// m_outStream << m_riState.IsShearPresent()   << zEndl;
    dst<<"0 0 0\n";

    for(const auto &poly : s.polymers){
        for(unsigned bead_id : poly.bead_ids){
            const auto &b=s.beads.at(bead_id);
            dst<<(poly.polymer_id+1); // Polymer ids are 1-based in Osprey
            dst<<" "<<poly.polymer_type; // Polymer types are 0-based
            dst<<" "<<(b.bead_id+1); // Bead ids are 1-based
            dst<<" "<<b.bead_type;
            dst<<" "<<s.bead_types[b.bead_type].r;
            dst<<" "<<b.x[0]<<" "<<b.x[1]<<" "<<b.x[2];
            dst<<" "<<b.x[0]<<" "<<b.x[1]<<" "<<b.x[2]; // TODO: what is unPBC XPos for us?
            dst<<" "<<b.v[0]<<" "<<b.v[1]<<" "<<b.v[2];
            dst<<" "<<b.f[0]<<" "<<b.f[1]<<" "<<b.f[2];
            dst<<"\n";
        }
    }

    if(!footer_strings.empty()){
        if(footer_strings[0]!="inclusive"){
            throw std::runtime_error("Invalid footer strings.");
        }
        for(const auto &s : footer_strings){
            dst<<s<<"\n";
        }
    }
}

struct osprey_bond_pattern_t
{
    std::string head;
    std::string tail;
    double kappa;
    double r0;
};
struct osprey_bond_pair_pattern_t
{
    std::string head;
    std::string mid;
    std::string tail;
    double kappa;
    double theta0;
};

struct osprey_patterns
{
    // Order is important
    std::vector<osprey_bond_pattern_t> bonds;
    std::vector<osprey_bond_pair_pattern_t> bond_pairs;

    std::unordered_map<std::pair<std::string,std::string>,osprey_bond_pattern_t> bonds_by_name;
    std::unordered_map<std::pair<std::string,std::string>,osprey_bond_pattern_t> bond_pairs_by_name;

    void add_bond(const std::string &head, const std::string &tail, double kappa, double r0)
    {
        auto it=bonds_by_name[{head,tail}];
        if(it!=bonds_by_name.end()){
            throw std::runtime_error("Repeated bond.");
        }
        it=bonds_by_name[{tail,head}];
        if(it!=bonds_by_name.end()){
            throw std::runtime_error("Repeated bond.");
        }
        bonds.push_back({head,tail,kappa,r0});
        bonds_by_name[{head,tail}]={head,tail,kappa,r0};
        bonds_by_name[{tail,head}]={head,tail,kappa,r0};
    }

    void add_bond_pair(const std::string &head, const std::string &tail, const std::string &tail, double kappa, double theta0)
    {
        auto it=bond_pairs_by_name[{head,mid,tail}];
        if(it!=bond_pairs_by_name.end()){
            throw std::runtime_error("Repeated bond pair.");
        }
        it=bond_pairs_by_name[{head,mid,tail}];
        if(it!=bond_pairs_by_name.end()){
            throw std::runtime_error("Repeated bond pair.");
        }
        bond_pairs.push_back({head,mid,tail,kappa,theta0});
        bond_pairs_by_name[{head,mid,tail}]={head,mid,tail,kappa,theta0};
        bond_pairs_by_name[{tail,mid,head}]={head,mid,tail,kappa,theta0};
    }
};

// The initial "dpd" line will already have been consumed
osprey_patterns parse_dmpci_to_osprey_patterns(std::istream &src)
{
    osprey_patterns patterns;

    // This is a dmpci file. Doesn't need to be fast
    std::string line;
    std::string keyword;
    while(std::getline(src, line)){
        src.str(line);
        if(! (src>>keyword) ){
            continue;
        }
        if(keyword == "Bond"){
            std::string head, tail;
            double kappa, r0;
            src >> head >> tail >> kappa >> r0;
            if(!src){
                throw std::runtime_error("Couldn parse Bond from dmpci file "+path+" : "+line);
            }
            
            patterns.add_bond(head, tail, kappa, r0);
        }else if(keyword == "BondPair"){
            std::string head, mid, tail;
            double kappa, theta0;
            src >> head >> mid >> tail >> kappa >> theta0;
            if(!src){
                throw std::runtime_error("Couldn parse Bond from dmpci file "+path+" : "+line);
            }
            
            patterns.add_bond_pair(head, mid, tail, theta0);
        }
    }
    return patterns;
}

// This assumes that the initial "inclusive" has already been consumed and parsed
osprey_patterns parse_restart_to_osprey_patterns(std::istream &src)
{
    // This can be slow, as the bulk of the bead data has already gone

    std::vector<std::string> bead_types;
    unsigned nBeadTypes;
    if(!(src >> nBeadTypes)){
        throw std::runtime_error("No bead type count in osprey template restart file "+path);
    }
    for(unsigned i=0; i<nBeadTypes; i++){
        unsigned bti;
        std::string name;
        src >> bti >> name;
        if(bti!=i){
            throw std::runtime_error("Incorrect bead type index in osprey template restart file "+path);
        }
        bead_types.push_back(name);
        
        double strength;
        for(unsigned j=0; j<2*nBeadTypes; j++){
            if(! (src >> strength) ){
                throw std::runtime_error("Error reading bead type interactions.");
            }
        }
    }

    /* This dance is because in restart states the bond and bond pair names
        are straight concatenations of bead types. So if we have
        bead types BA, AB, B, and AAB, then there is no way
        to tell what BAAB means. This is probably over-thinking, but
        short one-, two-, and three-letter bead-types are common.
    */
    std::unordered_map<std::string,std::vector<std::pair<std::string,std::string>>> bead_pairs;
    std::unordered_map<std::string,std::vector<std::tuple<std::string,std::string,std::string>>> bead_triples;
    for(unsigned i0=0; i0<nBeadTypes; i0++){
        auto n0=bead_types[i0];
        for(unsigned i1=0; i1<nBeadTypes; i1++){
            auto n1=bead_types[i1];

            auto n01=n0+n1;
            bead_pairs[n01].push_back({n0,n1});
            auto n10=n1+n0;
            if(n01 != n10){
                bead_pairs[n10].push_back({n1,n0});   
            }

            for(unsigned i2=0; i2<nBeadTypes; i2++){
                auto n2=bead_types[i2];

                auto n012=n0+n1+n2;
                bead_triples[n012].push_back({n0,n1,n2});
                auto n210=n2+n1+n0;
                if(n210 != n012){
                    bead_triples[n210].push_back({n2,n1,n0});
                }
            }
        }
    }

    unsigned nBonds;
    if(!(src >> nBonds)){
        throw std::runtime_error("No bond count in osprey template restart file.");
    }
    for(unsigned i=0; i<nBonds; i++){
        unsigned bi;
        std::string head, tail;
        double kappa, r0;
        src >> bi >> name >> kappa >> r0;
        if(bi!=i){
            throw std::runtime_error("Incorrect bead type index in osprey template restart file "+path);
        }
    }
}

std::vector<std::string> load_osprey_inclusive_restart_footer(
    const WorldState &state,
    const std::string &path
)
{
    osprey_patterns patterns;

    auto strip_trailing_ws=[](std::string &line) -> void
    {
        while(line.size()>0 && std::isspace(line.back())){
            line.pop_back();
        }
    };

    with_optional_gzip_istream(path, [&](std::istream &src){    

        std::string line;
        if(! std::getline(src, line)){
            throw std::runtime_error("Couldn't read first line from "+path);
        }

        strip_trailing_ws(line);
        if(line=="dpd"){
            patterns = parse_dmpci_to_osprey_patterns(src);
        }else{
            do{
                strip_trailing_ws(line);

                if(line=="inclusive"){
                    patterns = parse_restart_to_osprey_patterns(src);
                    break;
                }
            }while( std::getline( src , line ) );
            throw std::runtime_error("Couldn't treat "+path+" as either dmpci or osprey.");
        }
    });

    return patterns;
}

std::ostream &write_to_osprey_restart(
    std::ostream &dst,
    const WorldState &s,
    const std::string &dat_path={}
){
    std::vector<std::string> footer;

    if(!dat_path.empty()){
        footer=load_osprey_inclusive_restart_footer(dat_path);
    }

    return write_to_osprey_restart(dst, s, footer);
}

#endif
