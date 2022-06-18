#ifndef dpd_analysis_write_world_state_graph_as_graphml_hpp
#define dpd_analysis_write_world_state_graph_as_graphml_hpp

#include "dpd/analysis/walk_world_state_graph.hpp"

void write_world_state_graph_as_graphml(
    const WorldState &s,
    const std::function<bool(const Bead &)> &filter,
    std::ostream &dst
){
    dst<<
R"(<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"  
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="bt" for="node" attr.name="bead_type" attr.type="string" />
  <key id="pt" for="node" attr.name="polymer_type" attr.type="string" />
  <key id="po" for="node" attr.name="polymer_offset" attr.type="integer" />
  <key id="dist" for="edge" attr.name="dist" attr.type="double" />
  <key id="bnd" for="edge" attr.name="bonded" attr.type="boolean"><default>false</default></key>
  <key id="sp" for="edge" attr.name="same_poly" attr.type="boolean"><default>false</default></key>
  <graph id="G" edgedefault="undirected">
    
)";
    auto add_node=[&](const Bead &b)
    {
        dst<<"<node id='n"<<b.bead_id<<"'>";
        dst<<"<data key='bt'>"<<s.bead_types[b.bead_type].name<<"</data>";
        dst<<"<data key='pt'>"<<s.polymer_types[b.polymer_type].name<<"</data>";
        dst<<"<data key='po'>"<<b.polymer_offset<<"</data>";
        dst<<"</node>\n";
    };
    auto add_edge=[&](const Bead &a,const Bead &b, const vec3r_t &dx, double dr)
    {
        dst<<"<edge source='n"<<a.bead_id<<"' target='n"<<b.bead_id<<"'>";
        dst<<"<data key='dist'>"<<dr<<"</data>";
        if(a.polymer_id==b.polymer_id){
            bool is_bond=false;
            for(const Bond &bond : s.polymer_types[a.polymer_type].bonds){
                if( bond.bead_offset_head==a.polymer_offset && bond.bead_offset_tail==b.polymer_offset ){
                    is_bond=true;
                    break;
                }
                if( bond.bead_offset_head==b.polymer_offset && bond.bead_offset_tail==a.polymer_offset ){
                    is_bond=true;
                    break;
                }
            }
            if(is_bond){
                dst<<"<data key='bnd'>true</data>";
            }
            dst<<"<data key='sp'>true</data>";
        }
        dst<<"</edge>\n";
    };
    walk_world_state_graph(
        s, 2, filter, add_node, add_edge
    );

    dst<<
R"(
    </graph>
</graphml>
)";

}

#endif

