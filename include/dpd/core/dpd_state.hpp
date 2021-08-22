#ifndef dpd_state_hpp
#define dpd_state_hpp

#include <cstdint>
#include <vector>
#include <string>
#include <cassert>

#include "vec3.hpp"

// Here we generally assume that 32-bit is enough for indices and ids
// If Julian wants more than 4B beads, then... we need more money.

// This is a completely standardised hash-code, so that we get repeatable results across implementations.
/*  bbbb 1ppp pppp pppp  pppp pppp pppp pppp  : monomer, up to 2^27 instances
    bbbb 0ooo oooo pppp  pppp pppp pppp pppp  : polymer, up to 2^20 instances, and 127 beads per polymer
    b is the bead type index
*/
inline uint32_t bead_hash_construct(uint32_t bead_type, bool is_monomer, uint32_t polymer_id, uint32_t polymer_offset)
{
    assert( is_monomer ? (polymer_offset==0 && polymer_id<(1u<<27))
                                    : (polymer_offset<128 && polymer_id<(1u<<20)));
    uint32_t base=(uint32_t(polymer_offset)<<20) | polymer_id;
    base |= uint32_t(is_monomer)<<27;
    base |= bead_type << 28;
    return base;
}

inline bool bead_hash_is_monomer(uint32_t hash)
{ return hash&(1<<27); }

inline uint32_t bead_hash_get_bead_type(uint32_t hash)
{ return hash>>28; }

inline uint32_t bead_hash_get_polymer_offset(uint32_t hash)
{ return bead_hash_is_monomer(hash) ? 0 : ((hash>>20)&0x7F); }

inline uint32_t bead_hash_get_polymer_id(uint32_t hash)
{ return bead_hash_is_monomer(hash) ? (hash&0x7FFFFFF) : (hash&0xFFFFF); }

inline bool bead_hash_in_same_polymer(uint32_t h1, uint32_t h2)
{
    assert(h1!=h2);
    if( !(h1&h2&(1<<27)) ){
        return false;
    }
    return (h1&0xFFFFF) == (h2&0xFFFFF);
}

struct Bead
{
    // A unique bead id within the world. Bead ids must be contiguous, starting at zero.
    uint32_t bead_id = -1;

    // A unique polymer id within the world. Polymer ids must be contiguous, starting at zero.
    uint32_t polymer_id = -1;
    uint8_t polymer_offset = -1;

    // Technically these can be recovered from the polymer_id and polymer_offset,
    // but it is convenient and cache-friendly to have them here.
    uint8_t bead_type = -1;
    uint8_t polymer_type = -1;
    bool is_monomer;

    vec3r_t x;
    vec3r_t v;
    vec3r_t f;

    uint32_t get_bead_type() const
    { return bead_type; }

    // This is a completely standardised hash-code, so that we get repeatable results across implementations.
    /*  0000 1ppp pppp pppp  pppp pppp pppp pppp  : monomer, up to 2^27 instances
        0000 0ooo oooo pppp  pppp pppp pppp pppp  : polymer, up to 2^20 instances, and 127 beads per polymer
    */
    uint32_t get_hash_code() const
    {
        return bead_hash_construct(bead_type, is_monomer, polymer_id, polymer_offset);
    }
};

struct BeadType
{
    std::string name;
    double r;
    uint32_t id = -1;
};

struct Bond
{
    double kappa;
    double r0;
    uint32_t bead_offset_head = -1;
    uint32_t bead_offset_tail = -1;
};

struct BondPair
{
    double kappa;
    double theta0;
    uint32_t bond_offset_head = -1;
    uint32_t bond_offset_tail = -1;
};

struct PolymerType
{
    std::string name;
    std::vector<unsigned> bead_types;
    std::vector<Bond> bonds;
    std::vector<BondPair> bond_pairs;
    uint32_t polymer_id = -1;
};

struct Polymer
{
    std::vector<uint32_t> bead_ids;
    uint32_t polymer_id = -1;
    uint32_t polymer_type = -1;
};

struct InteractionStrength
{
    double conservative;
    double dissipative;
};

struct WorldState
{
    vec3r_t origin = {0, 0, 0}; // This is only for import/export. All coordinates are in [0,box)
    vec3r_t box = {4, 4, 4};  // Size of the box
    double lambda = 0.5;
    double t = 1.0;
    double dt = 0.05;
    uint64_t seed = 1;
    std::vector<InteractionStrength> interactions; // Array of bead_types.size()*bead_types.size(). Strictly symmetric
    std::vector<BeadType> bead_types;
    std::vector<PolymerType> polymer_types;

    std::vector<Polymer> polymers;
    std::vector<Bead> beads;
};

#endif