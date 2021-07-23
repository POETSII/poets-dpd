#ifndef bond_db_hpp
#define bond_db_hpp

#include "dpd_state.hpp"

#include <algorithm>

/*
    In this model we assume that there is a hookean bond between all
    beads, but things are very sparse. So for any random pair of
    beads we just return kappa of 0, and only for bonded beads do
    we return kappa>0.

    We also allow at most one angular bond per bead...
*/
struct bond_db_concept
{
    unsigned get_max_landing_slots() const;

    uint16_t get_key(uint16_t polymer_type, uint16_t polymer_offset) const;

    void get_bond_info(
        uint16_t home_key,
        float &kappa,     // kappa>1 for bond, kappa==0 for normal
        float &r0,        // If kappa==0 then this doesnt matter
        int &landing_zone // -1 if no landing zone
    ) const;

    bool get_bond_pair_info(
        uint16_t home_key,
        float &sin_theta0,
        float &sin_theta1,
        uint16_t &head_polymer_offset,
        uint16_t &tail_polymer_offset
    ) const;
};


struct bond_db_simple
{
public:
    //////////////////////////////////
    // Setup

    void clear()
    {
        m_max_landing_slots=0;
        m_nodes.clear();
        m_polymer_to_key.clear();
        m_bonds.clear();
    }

    void attach(const std::vector<PolymerType> &types)
    {
        clear();

        // Assign unique key to each polymer type member
        for(unsigned pti=0; pti<types.size(); pti++){
            const PolymerType &pt = types[pti];
            for(unsigned oi=0; oi<pt.bead_types.size(); oi++){
                m_polymer_to_key.push_back({pti,oi});
            }
        }
        
        // Loop through every bead in every polymer type, and allocate space
        for(unsigned pti=0; pti<types.size(); pti++){
            const PolymerType &pt = types[pti];
            for(unsigned oi=0; oi<pt.bead_types.size(); oi++){
                node_info node;
                node.bonds_begin=-1;
                node.bonds_end=-1;
                node.kappa=0;
                node.head_landing_zone=-1;
                node.tail_landing_zone=-1;
                node.head_polymer_offset=-1;
                node.tail_polymer_offset=-1;
                m_nodes.push_back(node);
            }
        }

        // Loop through every bond in every polymer
        std::vector<std::vector<bond_info>> node_to_bonds(m_nodes.size());
        for(unsigned pti=0; pti<types.size(); pti++){
            const PolymerType &pt = types[pti];
            for(unsigned bi=0; bi<pt.bonds.size(); bi++){
                const Bond &b = pt.bonds[bi];
                auto head_key = get_key(pt.polymer_id, b.bead_offset_head);
                auto tail_key = get_key(pt.polymer_id, b.bead_offset_tail);
                node_to_bonds[head_key].push_back({ tail_key, -1, b.kappa, b.r0 });
                node_to_bonds[tail_key].push_back({ head_key, -1, b.kappa, b.r0 });
            }
        }
        
    }



    unsigned get_max_landing_slots() const
    { return 2; } // We only support one angle per bead

    uint16_t get_key(uint16_t polymer_type, uint16_t polymer_offset) const
    {
        auto it=std::lower_bound(m_polymer_to_key.begin(), m_polymer_to_key.end(), {polymer_type, polymer_offset});
        assert(it->first==polymer_type && it->second==polymer_offset);
        return it-m_polymer_to_key.begin();
    }

    ////////////////////////////////////////////
    // Run-time

    void get_bond_info(
        uint16_t home_key,
        float &kappa,     // kappa>1 for bond, kappa==0 for normal
        float &r0,        // If kappa==0 then this doesnt matter
        int &landing_zone // -1 if no landing zone
    ) const 
    {
        const node_info &ni = m_nodes.at(home_key);
        const auto *begin=&m_bonds[ni.bonds_begin];
        const auto *end=&m_bonds[ni.bonds_end];

        while(begin!=end){
            if(begin->other==home_key){
                kappa=begin->kappa;
                r0=begin->r0;
                landing_zone=begin->landing_zone;
                return;
            }
            ++begin;
        }

        landing_zone=-1;
        kappa=0;
    }

    bool get_bond_pair_info(
        uint16_t home_key,
        float &kappa,
        float &sin_theta0,
        float &cos_theta0,
        uint16_t &head_landing_zone,
        uint16_t &tail_landing_zone,
        uint16_t &head_polymer_offset,
        uint16_t &tail_polymer_offset
    ) const
    {
        const node_info &ni = m_nodes.at(home_key);
        if(ni.kappa==0){
            return false;
        }
        kappa=ni.kappa;
        sin_theta0=ni.sin_theta0;
        cos_theta0=ni.cos_theta0;
        head_polymer_offset=ni.head_polymer_offset;
        tail_polymer_offset=ni.tail_polymer_offset;
        return true;
    }
private:
    struct bond_info
    {
        uint16_t other;
        int16_t landing_zone;
        float kappa;
        float r0;
    };

    struct node_info
    {
        uint16_t bonds_begin;
        uint16_t bonds_end;
        float kappa;
        float sin_theta0;
        float cos_theta0;
        uint16_t head_landing_zone;
        uint16_t tail_landing_zone;
        uint16_t head_polymer_offset;
        uint16_t tail_polymer_offset;
    };

    // Not really run-time things
    unsigned m_max_landing_slots;
    std::vector<std::pair<uint32_t,uint16_t>> m_polymer_to_key;

    // Run-time things
    std::vector<node_info> m_nodes;
    std::vector<bond_info> m_bonds;
};