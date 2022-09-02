#include "dpd/core/world_stream.hpp"

#include <iostream>

int main()
{
    double worst=0;

    vec3i_t delta;
    for(delta[0]=-10000; delta[0]<=+10000; delta[0]+=7){
        for(delta[1]=-10000; delta[1]<=+10000; delta[1]+=17){
            for(delta[2]=-10000; delta[2]<=+10000; delta[2]+=31){
                bead_delta_t enc{delta};

                vec3i_t got=enc.to_delta();
                double max_abs=0;
                for(int d=0; d<3; d++){
                    max_abs=std::max<double>(max_abs, std::abs(delta[d]));
                }
                double max_rel_err=0, max_abs_err=0;
                for(int d=0; d<3; d++){
                    max_abs_err=std::max<double>(max_abs_err, std::abs(delta[d]-got[d]));
                    max_rel_err=std::max<double>(max_rel_err, std::abs((delta[d]-got[d])/max_abs));
                }

                if(max_rel_err > worst){
                    std::cout<<delta<<" -> (s="<<enc.shift<<",dx="<<enc.dx<<","<<enc.dy<<","<<enc.dz<<") "<<got<<", max_abs="<<max_abs<<", max_rel="<<max_rel_err<<", max_abs="<<max_abs_err<<"\n";
                    worst=max_rel_err;
                    if(max_rel_err > 0.01){ // We have about 8 bit relative accuracy
                        std::cerr<<"Error budget exceeded.";
                        exit(1);
                    }
                }
            }
        }
    }
}