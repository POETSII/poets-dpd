#ifndef dpd_state_hpp
#define dpd_state_hpp

#include <cstdint>
#include <vector>
#include <string>
#include <cassert>
#include <unordered_map>

#include "vec3.hpp"

// Here we generally assume that 32-bit is enough for indices and ids
// If Julian wants more than 4B beads, then... we need more money.

// This is a completely standardised hash-code, so that we get repeatable results across implementations.
/*  bbbb 1ppp pppp pppp  pppp pppp pppp pppp  : monomer, up to 2^27 instances
    bbbb 0ooo ooop pppp  pppp pppp pppp pppp  : polymer, up to 2^21 instances, and 63 beads per polymer
    b is the bead type index
*/
struct BeadHash
{
    static const unsigned MONOMER_ID_BITS = 27;
    static const unsigned POLYMER_ID_BITS = 21;
    static const unsigned BEAD_TYPE_BITS = 32 - MONOMER_ID_BITS - 1;
    static const unsigned POLYMER_OFFSET_BITS = 32 - BEAD_TYPE_BITS - 1 - POLYMER_ID_BITS;

    uint32_t hash;

    BeadHash()
        : hash(0)
    {}

    explicit BeadHash(uint32_t raw)
        : hash(raw)
    {}

    static BeadHash construct(uint32_t bead_type, bool is_monomer, uint32_t polymer_id, uint32_t polymer_offset)
    {
        assert( is_monomer ? (polymer_offset==0 && polymer_id<(1u<<27))
                                        : (polymer_offset<64 && polymer_id<(1u<<21)));
        uint32_t base=(uint32_t(polymer_offset)<<21) | polymer_id;
        base |= uint32_t(is_monomer)<<27;
        base |= bead_type << 28;
        return BeadHash(base);
    }

    static BeadHash construct_checked(uint32_t bead_type, bool is_monomer, uint32_t polymer_id, uint32_t polymer_offset)
    {
        if(bead_type >= (1<<BEAD_TYPE_BITS)) throw std::runtime_error("bead_type is too large for hash.");
        if(is_monomer){
            if(polymer_id >= (1<<MONOMER_ID_BITS)) throw std::runtime_error("monomer id is too large for hash.");
        }else{
            if(polymer_offset >= (1<<POLYMER_OFFSET_BITS)) throw std::runtime_error("polymer offset is too large for hash.");
            if(polymer_id >= (1<<POLYMER_ID_BITS)) throw std::runtime_error("polymer id is too large for hash.");
        }
        return construct(bead_type, is_monomer, polymer_id, polymer_offset);
    }

    template<class T>
    static auto is_monomer(const T &hash)
    { return hash&(1<<27); }

    bool is_monomer() const
    { return is_monomer(hash); }

    template<class T>
    static auto get_bead_type(const T &hash)
    { return hash>>28; }

    uint32_t get_bead_type() const
    { return get_bead_type(hash); }

    inline uint32_t get_bead_type() const
    { return hash>>28; }

    template<class T>
    static auto get_bead_type(const T &hash)
    { return hash>>28; }

    inline uint32_t get_bead_type() const
    { return get_bead_type(hash); }


    inline uint32_t get_polymer_offset() const
    { return is_monomer() ? 0 : ((hash>>21)&0x3F); }

    inline uint32_t get_polymer_id() const
    { return is_monomer() ? (hash&0x7FFFFFF) : (hash&0x1FFFFF); }

    inline bool in_same_polymer(const BeadHash &h2) const
    {
        assert(hash!=h2.hash);
        if( !(hash&h2.hash&(1<<27)) ){
            return false;
        }
        return (hash&0x1FFFFF) == (h2.hash&0x1FFFFF);
    }

    inline BeadHash reduced_hash() const
    { return BeadHash{hash&uint32_t(0xFFFFFFFul)}; }

    // This cannot work out what the bead type is, so will return the 0 for the bead type
    inline BeadHash make_reduced_hash_from_polymer_offset(unsigned polymer_offset) const
    {
        assert(!is_monomer());
        assert(polymer_offset < 128);
        return BeadHash{ (hash & 0x1FFFFF) | (polymer_offset<<21) }; 
    }

    inline bool reduced_equals(const BeadHash &h2) const
    {
        return (hash&0x0FFFFFFFul)==(h2.hash&0x0FFFFFFFul);
    }

    inline bool operator==(const BeadHash &o) const
    { return hash==o.hash; }

    inline bool operator!=(const BeadHash &o) const
    { return hash!=o.hash; }
};

namespace std
{
    template<>
    struct hash<BeadHash>
    {
        size_t operator()(BeadHash h) const
        {
            return h.hash;
        }
    };
};

static const float MIN_DISTANCE_CUTOFF = 0.000000001;
static const float MIN_DISTANCE_CUTOFF_SQR = MIN_DISTANCE_CUTOFF * MIN_DISTANCE_CUTOFF;

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

    BeadHash get_hash_code() const
    {
        return BeadHash::construct(bead_type, is_monomer, polymer_id, polymer_offset);
    }

    BeadHash get_hash_code_checked() const
    {
        return BeadHash::construct(bead_type, is_monomer, polymer_id, polymer_offset);
    }

    void set_hash_code(BeadHash h)
    {
        bead_type=h.get_bead_type();
        is_monomer=h.is_monomer();
        polymer_id=h.get_polymer_id();
        polymer_offset=h.get_polymer_offset();
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

    bool operator==(const Bond &o) const
    { return kappa==o.kappa && r0==o.r0 && bead_offset_head==o.bead_offset_head && bead_offset_tail==o.bead_offset_tail; }
};

struct BondPair
{
    double kappa;
    double theta0;
    uint32_t bond_offset_head = -1;
    uint32_t bond_offset_tail = -1;

    bool operator==(const BondPair &o) const
    { return kappa==o.kappa && theta0==o.theta0 && bond_offset_head==o.bond_offset_head && bond_offset_tail==o.bond_offset_tail; }
};

struct PolymerType
{
    std::string name;
    std::vector<unsigned> bead_types;
    std::vector<Bond> bonds;
    std::vector<BondPair> bond_pairs;
    uint32_t polymer_id = -1;

    bool operator==(const PolymerType &o) const
    {
        return name==o.name && bead_types==o.bead_types && bonds==o.bonds && bond_pairs==o.bond_pairs && polymer_id==o.polymer_id; 
    }
};

struct Polymer
{
    std::vector<uint32_t> bead_ids;
    uint32_t polymer_id = -1;
    uint32_t polymer_type = -1;

    bool operator==(const Polymer &o) const
    { return bead_ids==o.bead_ids && polymer_type==o.polymer_type && polymer_id==o.polymer_id; }
};

struct InteractionStrength
{
    double conservative;
    double dissipative;

    bool operator==(const InteractionStrength &o) const
    { return conservative==o.conservative && dissipative==o.dissipative; }
};

struct WorldState
{
    vec3r_t origin = {0, 0, 0}; // This is only for import/export. All coordinates are in [0,box)
    vec3r_t box = {4, 4, 4};  // Size of the box
    double lambda = 0.5;
    long t =0;
    double dt = 0.05;
    uint64_t seed = 1;
    std::vector<InteractionStrength> interactions; // Array of bead_types.size()*bead_types.size(). Strictly symmetric
    std::vector<BeadType> bead_types;
    std::vector<PolymerType> polymer_types;

    std::vector<Polymer> polymers;
    std::vector<Bead> beads;

    uint32_t bead_hash_to_id(const BeadHash &hash)
    {
        return polymers.at(hash.get_polymer_id()).bead_ids.at(hash.get_polymer_offset());
    }
};

#endif
