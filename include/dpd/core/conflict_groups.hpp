#ifndef make_conflict_groups_hpp
#define make_conflict_groups_hpp

#include <unordered_set>
#include "dpd/core/vec3.hpp"

// Turns locations into clusters of points, and assigns those clusters to
// rounds. All clusters within a round should be seperated by at least
// seperation in all dimensions (including wrap-around). The clusters
// are also expected to be compact, to maximise locality.
// The result is a vector of vector of vector, nested as:
//  - Groups in conflict group
//  - Clusters in Group
//  - Points in cluster
auto make_conflict_groups(vec3i_t dims, int seperation) -> std::vector<std::vector<std::vector<vec3i_t>>>
{
    if( (dims[0]%4) || (dims[1]%4) || (dims[2]%4)){
        throw std::runtime_error("make_conflict_groups : Currently dims must be a multiple of 4 in each direction.");
    }
    if(seperation>2){
        throw std::runtime_error("Don't know how to do seperation>2");
    }

    std::unordered_set<unsigned> seen;
    std::vector<std::vector<std::vector<vec3i_t>>> conflict_groups;

    for(unsigned gx=0; gx<4; gx+=2){
        for(unsigned gy=0; gy<4; gy+=2){
            for(unsigned gz=0; gz<4; gz+=2){

                // This gives the origin of a 2x2x2 cube within a 4x4x4 block
                std::vector<std::vector<vec3i_t> > group;
                // Loop over all super cells within group
                for(unsigned ix=gx; ix<(unsigned)m_dims[0]; ix+=4){
                    for(unsigned iy=gy; iy<(unsigned)m_dims[1]; iy+=4){
                        for(unsigned iz=gz; iz<(unsigned)m_dims[2]; iz+=4){
                            
                            std::vector<vec3i_t> cluster;
                            for(int lx=0; lx<2; lx++){
                                for(int ly=0; ly<2; ly++){
                                    for(int lz=0; lz<2; lz++){
                                        vec3i_t pos{int(ix+lx),int(iy+ly),int(iz+lz)});
                                        if(!seen.insert(pos).second){
                                            throw std::runtime_error("Duplicate");
                                        }
                                        cluster.push_back(pos);
                                    }
                                }
                            }
                            group.push_back(std::move(sc));
                        }
                    }
                }
                conflict_groups.push_back(std::move(group));
            }
        }
    }

    return conflict_groups;
}

#endif