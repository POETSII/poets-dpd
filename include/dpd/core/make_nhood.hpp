#ifndef make_nhood_hpp
#define make_nhood_hpp

#include "dpd/core/vec3.hpp"
#include <vector>
#include <unordered_set>
#include <unordered_map>

enum struct NHoodType
{
    Full =      1,
    Forwards =  2,
    Backwards = 3
};

template<class TFunc>
inline void for_each_point_in_box(const vec3i_t &begin, const vec3i_t &end, const TFunc &func)
{
    vec3i_t x;
    for(x[0]=begin[0]; x[0] < end[0] ;x[0]++){
        for(x[1]=begin[1]; x[1] < end[1] ;x[1]++){
            for(x[2]=begin[2]; x[2] < end[2] ;x[2]++){
                func(x);
            }
        }
    }
}

inline std::vector<vec3i_t> make_relative_nhood_full(bool exclude_centre=false)
{
    std::vector<vec3i_t> res;
    for_each_point_in_box({-1,-1,-1},{+2,+2,+2},[&](const vec3i_t &dir){
        if(!exclude_centre || dir[0]!=0 || dir[1]!=0 || dir[2]!=0){
            res.push_back(dir);
        }
    });
    return res;
}

inline std::vector<vec3i_t> make_relative_nhood_forwards(bool exclude_centre=false)
{
    std::unordered_set<vec3i_t> res;
    for(auto dir : make_relative_nhood_full(exclude_centre)){
        if(res.find(-dir)==res.end()){
            res.insert(dir);
        }
    }
    return {res.begin(),res.end()};
}

inline std::vector<vec3i_t> make_relative_nhood_inverse(const std::vector<vec3i_t> &fwds)
{
    std::vector<vec3i_t> res{fwds};
    for(auto &dir : res){
        dir = -dir;
    }
    return res;
}

inline std::vector<vec3i_t> make_relative_nhood(NHoodType nhood, bool exclude_centre)
{
    switch(nhood){
        case NHoodType::Full: return make_relative_nhood_full(exclude_centre);
        case NHoodType::Forwards: return make_relative_nhood_forwards(exclude_centre);
        case NHoodType::Backwards: return make_relative_nhood_inverse(make_relative_nhood_forwards(exclude_centre));
        default: throw std::logic_error("Invalid nhood type.");
    }
}

inline std::vector<vec3i_t> make_absolute_nhood(const std::vector<vec3i_t> &rel_nhood, vec3i_t size, vec3i_t centre)
{
    std::vector<vec3i_t> res;
    res.reserve(rel_nhood.size());
    for(auto dir : rel_nhood){
        vec3i_t tmp=centre+dir;
        for(int i=0; i<3; i++){
            tmp[i] = (tmp[i]+size[i])%size[i];
        }
        res.push_back(tmp);
    }
    return res;
}

inline std::unordered_map<vec3i_t,std::vector<vec3i_t>> make_absolute_nhood_map(const std::vector<vec3i_t> &rel_nhood, vec3i_t size)
{
    std::unordered_map<vec3i_t,std::vector<vec3i_t>> res;
    res.reserve(size[0]*size[1]*size[2]);
    for_each_point_in_box({0,0,0}, size, [&](const vec3i_t &centre){
        res[centre]=std::move(make_absolute_nhood(rel_nhood, size, centre));
    });
    return res;
}

#endif