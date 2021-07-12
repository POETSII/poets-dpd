#ifndef dpd_state_builder_hpp
#define dpd_state_builder_hpp

#include "dpd_state.hpp"

#include "vec3.hpp"

#include <variant>

struct WorldStateBuilderError
    : std::runtime_error
{
    WorldStateBuilderError(const std::string &s)
        : std::runtime_error(s)
    {}
};

class WorldStateBuilder
{
private:
    double m_default_interaction_conservative = 1;
    double m_default_interaction_dissipative = 0.1;
    WorldState m_world;

    void require(bool cond, const std::string &msg)
    {
        if(!cond){
            throw WorldStateBuilderError(msg);
        }
    }
public:
    WorldStateBuilder(vec3r_t dims)
    {
        for(unsigned i=0; i<3; i++){
            require(round(dims[i])==dims[i], "Dimensions must be integer.");
            m_world.box[i]=dims[i];
        }
    }

    WorldState &data()
    { return m_world; }

    WorldState extract()
    { return std::move(m_world); }

    using index_or_name_t = std::variant<unsigned,std::string>;

    unsigned add_bead_type(std::string name)
    {
        for(const auto &b : m_world.bead_types){
            require(b.name!=name, "Bead type name already exists.");
        }
        unsigned index=m_world.bead_types.size();
        m_world.bead_types.push_back({
            name, 0.5, index
        });

        unsigned n=m_world.bead_types.size();
        std::vector<InteractionStrength> newInteractions(n*n);
        for(unsigned i=0; i<n; i++){
            for(unsigned j=0; j<n; j++){
                if(i<(n-1) && j<(n-1)){
                    newInteractions[i*n+j] = m_world.interactions[i*(n-1)+j];
                }else{
                    newInteractions[i*n+j] = {m_default_interaction_conservative,m_default_interaction_dissipative};
                }
            }
        }
        std::swap(newInteractions, m_world.interactions);

        return index;
    }

    void set_interaction_strength(const index_or_name_t &a, const index_or_name_t &b, double conservative, double dissipative)
    {
        unsigned ia=get_bead_type(a);
        unsigned ib=get_bead_type(b);
        m_world.interactions.at(ia*m_world.bead_types.size()+ib)={conservative,dissipative};
        m_world.interactions.at(ib*m_world.bead_types.size()+ia)={conservative,dissipative};
    }

    unsigned get_bead_type(const index_or_name_t &n)
    {
        if(n.index()==0){
            unsigned index=std::get<0>(n);
            require(index < m_world.bead_types.size(), "Invalid bead type index.");
            return index;
        }else{
            for(const auto &b : m_world.bead_types){
                if(b.name==std::get<1>(n)){
                    return b.id;
                }
            }
            require(false, "Invalid bead type name.");
            return -1; // Not reachable
        }
    }

    unsigned get_polymer_type(const index_or_name_t &n)
    {
        if(n.index()==0){
            unsigned index=std::get<0>(n);
            require(index < m_world.polymer_types.size(), "Invalid polymer type index.");
            return index;
        }else{
            for(const auto &b : m_world.polymer_types){
                if(b.name==std::get<1>(n)){
                    return b.polymer_id;
                }
            }
            require(false, "Invalid polymer type name.");
            return -1; // Not reachable
        }
    }

    unsigned add_polymer_type(
        std::string name,
        const std::vector<index_or_name_t> &bead_types,
        const std::vector<Bond> &bonds,
        const std::vector<BondPair> &bond_pairs
    ){
        for(const auto &p : m_world.polymer_types){
            require(p.name!=name, "Polymer type name already exists.");
        }
        unsigned index=m_world.polymer_types.size();
        PolymerType pt;
        pt.name=name;
        pt.polymer_id=index;
        for(const auto &bt : bead_types){
            pt.bead_types.push_back(get_bead_type(bt));
        }
        for(const auto &b : bonds){
            require(b.bead_offset_head<bead_types.size(), "Invalid bead index.");
            require(b.bead_offset_tail<bead_types.size(), "Invalid bead index.");
            pt.bonds.push_back(b);
        }
        for(const auto &bp : bond_pairs){
            require(bp.bond_offset_head<pt.bonds.size(), "Invalid bond index.");
            require(bp.bond_offset_tail<pt.bonds.size(), "Invalid bond index.");
            pt.bond_pairs.push_back(bp);
        }
        m_world.polymer_types.push_back(pt);
        return index;
    }

    struct bead_info_t
    {
        vec3r_t x;
        vec3r_t v{0,0,0};
        vec3r_t f{0,0,0};
    };

    unsigned add_polymer(index_or_name_t poly_type,
        const std::vector<bead_info_t> &beads,
        bool wrap=false
    ){
        unsigned poly_type_index=get_polymer_type(poly_type);
        const auto &pt=m_world.polymer_types.at(poly_type_index);
        require(pt.bead_types.size() == beads.size(), "Wrong number of beads for polymer type.");

        unsigned polymer_id=m_world.polymers.size();

        std::vector<uint32_t> bead_ids;
        for(unsigned i=0; i<beads.size(); i++){
            vec3r_t x=beads[i].x;
            if(wrap){
                for(unsigned j=0; j<3; j++){
                    while(x[i] < 0){ x[i] += m_world.box[j]; }
                    while(x[i] >= m_world.box[j] ){ x[i] -= m_world.box[j]; }
                }
            }else{
                for(unsigned j=0; j<3; j++){
                    require(0 <= beads[i].x[j], "Bead out of box.");
                    require(m_world.box[j] > beads[i].x[j], "Bead out of box.");
                }
            }

            Bead b;
            b.bead_id=m_world.beads.size();
            bead_ids.push_back(b.bead_id);
            b.bead_type=pt.bead_types[i];
            b.polymer_type=poly_type_index;
            b.polymer_id=polymer_id;
            b.polymer_offset=i;
            b.x=x;
            b.v=beads[i].v;
            b.f=beads[i].f;
            m_world.beads.push_back(b);
        }

        m_world.polymers.push_back({std::move(bead_ids), polymer_id, poly_type_index});
        return polymer_id;
    }

    unsigned add_monomer(index_or_name_t poly_type, const vec3r_t &x, const vec3r_t &v=vec3r_t(), const vec3r_t &f=vec3r_t())
    {
        return add_polymer(poly_type, {{x,v,f}});
    }

    
};

#endif