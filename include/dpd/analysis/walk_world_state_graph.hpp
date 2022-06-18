#ifndef dpd_core_enumerate_connections_hpp
#define dpd_core_enumerate_connections_hpp

#include "dpd/core/dpd_state.hpp"
#include "dpd/core/make_nhood.hpp"

#include <functional>
#include <vector>

void walk_world_state_graph(
    const WorldState &state,
    int max_r,
    const std::function<bool(const Bead &a)> &filter_out,
    const std::function<void(const Bead &a)> &add_vertex,
    const std::function<void(const Bead &a, const Bead &b, const vec3r_t &dx, double dr )> &add_edge
){
    assert(max_r>0);

    vec3i_t dims{state.box};
    for(int d=0; d<3; d++){
        if( (dims[d]%max_r) ){
            throw std::runtime_error("enumerate_connections : max_r must divide state.dims");
        }
        if( dims[d] / max_r <= 3){
            // This is needed for wrap-around checking
            throw std::runtime_error("enumerate_connections : max_r must divide state.dims at least 4 times");
        }
        dims[d] /= max_r;
    }

    std::vector<std::vector<const Bead *>> cells;
    cells.resize(dims[0]*dims[1]*dims[2]);


    auto pos_to_cell_index=[&](vec3i_t ix)
    {
        for(int d=0; d<3; d++){
            if(ix[d]<0){
                ix[d]+=dims[d];
            }else if(ix[d]>=dims[d]){
                ix[d]-=dims[d];
            }
        }
        return ix[0] + ix[1]*dims[0] + ix[2]*dims[0]*dims[1];
    };

    auto cell_index_to_pos=[&](unsigned index) -> vec3i_t
    {
        unsigned oindex=index;
        vec3i_t res;
        res[0]=index%dims[0];
        index/=dims[0];
        res[1]=index%dims[1];
        index/=dims[2];
        res[2]=index;
        assert(oindex==pos_to_cell_index(res));
        return res;
    };

    auto bead_to_cell_index=[&](const Bead &b)
    {
        vec3i_t ix=vec3_floor(b.x) / max_r;
        auto index= pos_to_cell_index(ix);
        std::cerr<<"index "<<index<<" -> ix="<<ix<<"\n";
        return index;
    };


    double wrap_tol=3*max_r;

    auto interact=[&](const Bead &a, const Bead &b)
    {
        vec3r_t dx=a.x-b.x;
        for(int d=0; d<3; d++){
            if(dx[d] < -wrap_tol){
                dx[d] += state.box[d];
            }
            if(dx[d] > wrap_tol){
                dx[d] -= state.box[d];
            }
        }
        double r=dx.l2_norm();
        if(r < max_r){
            add_edge(a, b, dx, r);
        }
    };

    for(const auto &b : state.beads){
        if(filter_out(b)){
            continue;
        }
        unsigned index=bead_to_cell_index(b);
        cells.at(index).push_back(&b);
        add_vertex(b);
    }

    auto rel_nhood=make_relative_nhood_forwards(true);
    for(unsigned home_index=0; home_index<cells.size(); home_index++){
        const auto &home_cell=cells[home_index];
        auto home_pos=cell_index_to_pos(home_index);
        for(int i=0; i+1<home_cell.size(); i++){
            for(int j=i+1; j<home_cell.size(); j++){
                interact(*home_cell[i], *home_cell[j]);
            }
        }
        for(auto dir : rel_nhood){
            auto other_index=pos_to_cell_index( home_pos + dir );
            const auto &other_cell=cells[other_index];
            for(int i=0; i<home_cell.size(); i++){
                for(int j=i+1; j<other_cell.size(); j++){
                    interact(*home_cell[i], *other_cell[j]);
                }
            }
        }
    }
}

#endif
