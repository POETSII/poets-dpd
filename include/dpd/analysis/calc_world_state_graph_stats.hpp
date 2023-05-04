#error "Work in progress"

#include <vector>

#include "dpd/core/dpd_state.hpp"

#include "dpd/core/json_helper.hpp"

#include "robin_hood.h"

/*
We collect the following observables, averaged over K timesteps:
- bead_class:
    - Histogram bead_class * 32. Number of beads of bead_class within [0,1/16), [1/16,2/16), ... [ 31/16, 32/16)
    - average squared speed (temperature)
- average squared speed (temperature)

Bead class is a linear index identifying position within polymer 

*/
struct world_state_graph_stats
{
    robin_hood::unordered_flat_map<uint32_t,uint32_t> polymer_id_mul256_plus_offset_to_class;
    robin_hood::unordered_flat_map<uint32_t,uint32_t> bead_hash_to_class;

    world_state_graph_stats(const WorldState &s)
    {
        unsigned offset=0;
        for(const PolymerType &pt : s.polymer_types){
            if(pt.bead_types.size() >= 256){
                throw std::runtime_error("Polymers are too long/big.");
            }
            for(unsigned o=0; i<pt.bead_types.size(); i++){
                polymer_id_mul256_plus_offset_to_class[pt.polymer_id*256+o] = offset++;
        }

        bead_hash_to_class.reserve(s.beads.size());
        
        for(const Bead &b : s.beads){
            BeadHash h=b.get_hash_code();
            unsigned polymer_type=s.polymers.at(h.get_polymer_id()).polymer_type;
            unsigned polymer_offset=h.get_polymer_offset();
            bead_hash_to_class[h.hash] = polymer_id_mul256_plus_offset_to_class[ polymer_type*256 + polymer_offset ];
        }
    }

    unsigned get_bead_class(const BeadHash &h)
    {
        return bead_hash_to_class.at(h.hash);
    }

    unsigned get_quant_bead_radius(double dr)
    {
        assert(0<dr && dr <2);
        return floor(dr * 16);
    }

    unsigned get_bead_class_quant_radius(const BeadHash &h, double dr)
    {
        return get_bead_class(h) * get_quant_bead_radius(dr);
    }

    struct per_bead_class_stats
    {
        unsigned index; // Bead type index

        unsigned n;      // Number of beads of this type
        vec3r_t sum_v2;  // Sum of velocity squared for all beads of this type
    
        // Histogram of adjacency by bead_type then number of beads of that type
        std::vector<uint64_t> adjacency;
    };

    struct timeslice_stats
    {
        unsigned t;
        unsigned K;
        std::vector<per_bead_class_stats> bead_class_stats;
    };

    void on_bead(const Bead &a)
    {
        auto &bs = bead_class_stats.at(get_bead_class(a.get_hash_code()));
        bs.n += 1;
        bs.sum_v2 += a.v.dot_self();
    }

    void on_dpd_bond(const Bead &a, const Bead &b, const vec3r_t &dx, double dr )
    {
        auto &bs = bead_class_stats.at(get_bead_class(a.get_hash_code()));
        bs.adjacency.at(get_bead_class_quant_radius(b, dr));
    }



};