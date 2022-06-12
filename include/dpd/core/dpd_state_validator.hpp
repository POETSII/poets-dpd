#ifndef dpd_state_validator_hpp
#define dpd_state_validator_hpp

#include "dpd_state.hpp"

#include "dpd/core/vec3.hpp"

#include <unordered_set>
#include <unordered_map>
#include <iostream>

class ValidationFailedError
    : public std::runtime_error
{
public:
    ValidationFailedError(const std::string &s)
        : std::runtime_error(s)
    {}
};

inline void validate(const WorldState &s, double max_r=1, double max_dist_per_step=0.5)
{
    auto require=[](bool cond, const char *msg)
    {
        if(!cond){
            throw ValidationFailedError(msg);
        }
    };
#define REQUIRE(cond) require((cond), "WorldState validation failed : " #cond)

    REQUIRE(0 <= s.lambda && s.lambda <= 1);
    REQUIRE(0<s.dt && s.dt < 1);
    REQUIRE(s.t >= 0);

    for(unsigned i=0; i<3; i++){
        REQUIRE( s.box[i] == round(s.box[i] ) );
        REQUIRE( s.box[i] >= 4 );
    }
    
    std::unordered_set<std::string> bead_type_names;
    for(const BeadType &bt : s.bead_types){
        REQUIRE(bt.id==bead_type_names.size());
        REQUIRE(  bead_type_names.insert(bt.name).second );
        REQUIRE( 0 < bt.r && bt.r <= 0.5 );
    }

    REQUIRE(s.interactions.size() == s.bead_types.size() * s.bead_types.size());
    for(unsigned i=0; i<s.interactions.size(); i++){
        unsigned a=i / s.bead_types.size();
        unsigned b=i % s.bead_types.size();
        unsigned j = s.bead_types.size() * b + a;
        REQUIRE(s.interactions.at(i).conservative == s.interactions.at(j).conservative);
        REQUIRE(s.interactions.at(i).dissipative == s.interactions.at(j).dissipative);
        REQUIRE(0<= s.interactions.at(i).conservative);
        REQUIRE(0<= s.interactions.at(i).dissipative);
    }

    std::unordered_set<std::string> polymer_type_names;
    for(const PolymerType &pt : s.polymer_types){
        REQUIRE(pt.polymer_id==polymer_type_names.size());
        REQUIRE(  polymer_type_names.insert(pt.name).second );
        for(auto bt : pt.bead_types){
            REQUIRE(bt < s.bead_types.size());
        }
        for(const auto &b : pt.bonds){
            REQUIRE(b.bead_offset_head < pt.bead_types.size());
            REQUIRE(b.bead_offset_tail < pt.bead_types.size());
            REQUIRE(b.bead_offset_head!=b.bead_offset_tail);
            REQUIRE(0 < b.kappa);
            REQUIRE(0 < b.r0 && b.r0 < 1);
        }
        for(const auto &bp : pt.bond_pairs){
            REQUIRE(bp.bond_offset_head < pt.bonds.size());
            REQUIRE(bp.bond_offset_tail < pt.bonds.size());
            REQUIRE(bp.bond_offset_head != bp.bond_offset_tail);
            REQUIRE(pt.bonds.at(bp.bond_offset_head).bead_offset_tail == pt.bonds.at(bp.bond_offset_tail).bead_offset_head);
            REQUIRE(0 < bp.kappa);
            REQUIRE(0 <= bp.theta0 && bp.theta0 <= 3.1416);
        }
    }

    std::vector<bool> polymers_seen;
    polymers_seen.resize(s.polymers.size(), false);
    std::vector<bool> beads_seen;
    beads_seen.resize(s.beads.size(), false);
    for(const auto & p : s.polymers){
        REQUIRE(!polymers_seen.at(p.polymer_id));
        polymers_seen[p.polymer_id]=true;
        REQUIRE(p.polymer_type < s.polymer_types.size());
        const auto &pt = s.polymer_types.at(p.polymer_type);
        REQUIRE(p.bead_ids.size() == pt.bead_types.size());
        for(unsigned i=0; i<p.bead_ids.size(); i++){
            unsigned bid=p.bead_ids[i];
            const auto &b = s.beads.at(bid);
            REQUIRE( b.bead_type == pt.bead_types[i] );
            REQUIRE( (p.bead_ids.size()==1) == b.is_monomer);
            REQUIRE(b.polymer_offset==i);
            REQUIRE( !beads_seen.at(bid));
            beads_seen[bid]=true;
        }

        for(const auto &b : pt.bonds){
            const auto &h = s.beads.at( p.bead_ids.at( b.bead_offset_head) );
            const auto &t = s.beads.at( p.bead_ids.at( b.bead_offset_tail) );
            vec3r_t distance=h.x - t.x;
            for(int i=0; i<3; i++){
                if(std::abs(distance[i]+s.box[i]) < std::abs(distance[i])){
                    distance[i] += s.box[i];
                }
                if(std::abs(distance[i]-s.box[i]) < std::abs(distance[i])){
                    distance[i] -= s.box[i];
                }
            }
            double r=distance.l2_norm();
            REQUIRE( r < max_r);
        }
    }

    std::unordered_map<BeadHash,const Bead *> seen_hashes;

    // Max speed allowed is half a box per time-step
    double max_v = max_dist_per_step / s.dt;
    // Max force would cause speed to increase by max_v in one time-step
    double max_f = max_v / s.dt;
    // or distance to change by more than 0.5 in one time-step
    max_f = std::min(max_f, 2*max_dist_per_step/(s.dt*s.dt));
    // TODO : double-check those velocity and force thresholds.
    for(const auto & b : s.beads){
        for(unsigned i=0; i<3; i++){
            REQUIRE( 0 <= b.x[i] );
            double bx=b.x[i], sb=s.box[i];
            if(bx >= sb){
                std::cerr<<"  bead "<<i<<", x="<<b.x<<", dims="<<s.box<<"\n";
            }
            REQUIRE(b.x[i] < s.box[i]);

            if(!(-max_v <= b.v[i] && b.v[i] <= max_v)){
                std::cerr<<"  bead "<<i<<", v="<<b.v<<", max_v="<<max_v<<"\n";
                exit(1);
            }
            REQUIRE( -max_v <= b.v[i] && b.v[i] <= max_v);

            if(!(-max_f <= b.f[i] && b.f[i] <= max_f)){
                std::cerr<<"  bead "<<i<<", f="<<b.f<<", max_f="<<max_f<<"\n";
                exit(1);
            }
            REQUIRE( -max_f <= b.f[i] && b.f[i] <= max_f );
        }

        auto it=seen_hashes.find(b.get_hash_code());
        if(it!=seen_hashes.end()){
            std::cerr<<"  hash collision : "<<it->second->bead_id<<", "<<b.bead_id<<"\n";
        }
        REQUIRE( seen_hashes.insert({b.get_hash_code(),&b}).second);
    }


#undef REQUIRE
}

#endif