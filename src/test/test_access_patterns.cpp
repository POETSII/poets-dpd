#include "dpd/core/make_nhood.hpp"

#include <unordered_set>
#include <tuple>
#include <iostream>

struct stats
{
    unsigned total=0;
    unsigned misses=0;
    std::vector<unsigned> access_distance;
};

stats count_cache(
    const std::vector<vec3i_t> &rel_nhood,
    const std::vector<vec3i_t> &order
){
    std::unordered_map<vec3i_t,unsigned> seen;
    stats res;

    unsigned t=0;
    for(auto base : order){
        for(auto dir : rel_nhood){
            auto it=seen.find(base+dir);
            if(it==seen.end()){
                res.misses++;
                seen[base+dir]=t;
            }else{
                unsigned distance=t-it->second;
                if(distance >= res.access_distance.size()){
                    res.access_distance.resize(distance+1, 0);
                }
                res.access_distance[distance]++;
                it->second=t;
            }
            res.total++;
            t++;
        }
    }
    return res;
}

std::vector<vec3i_t> make_order(int d0, int w0, int d1, int w1, int d2, int w2)
{
    vec3i_t point;
    std::vector<vec3i_t> res;
    for(int i0=0; i0<w0; i0++){
        point[d0]=i0;
        for(int i1=0; i1<w1; i1++){
            point[d1]=i1;
            for(int i2=0; i2<w2; i2++){
                point[d2]=i2;
                res.push_back(point);
            }
        }
    }
    return res;
}

void test(const std::vector<vec3i_t> &nhood, int d0, int w0, int d1, int w1, int d2, int w2)
{
    auto order=make_order(d0, w0, d1, w1, d2, w2);
    auto res=count_cache(nhood, order);
    std::cout<<d0<<","<<w0<<", "<<d1<<","<<w1<<", "<<d2<<","<<w2<<", total="<<res.total<<", rate="<<res.misses/(double)res.total<<",   ";
    double sum=0;
    for(unsigned i=0; i<res.access_distance.size(); i++){
        sum += i * res.access_distance[i];
    }
    std::cout<<"mean="<<sum/res.total<<", max="<<res.access_distance.size()-1<<"\n";
    std::cout<<"\n";
}

int main()
{
    auto rel=make_relative_nhood_forwards(false);
    for(auto p : rel){
        std::cerr<<p<<"\n";
    }

    test(rel, 0,16, 1,2, 2,2);
    test(rel, 0,16, 2,2, 1,2);

    test(rel, 1,16, 0,2, 2,2);
    test(rel, 1,16, 2,2, 0,2);

    test(rel, 2,16, 0,2, 1,2);
    test(rel, 2,16, 1,2, 0,2);
}