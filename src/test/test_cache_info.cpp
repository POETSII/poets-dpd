#include "dpd/core/cache_info.hpp"

#include <iostream>

int main()
{
    auto cpu0=cache_info::build_cpu_cache_info(0);
    for(unsigned i=0; i<cpu0.size(); i++){
        std::cout<<cpu0[i].level<<", "<<cpu0[i].size<<"\n";
    }

    auto sys=cache_info::build_cpu_sets();
    for(auto &lev : sys.levels){
        std::cerr<<lev.level<<", node_size="<<lev.node_size<<", total_size="<<lev.total_size<<" :";
        for(auto &c : lev.nodes){
            std::cerr<<" "<<c;
        }
        std::cerr<<"\n";
    }
}
