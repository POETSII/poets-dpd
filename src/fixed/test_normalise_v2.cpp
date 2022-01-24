#include "ap_fixed_vec3.hpp"

#include <random>
#include <iostream>

int main()
{
    std::mt19937_64 urng;

    double max_err=0;
    double max_cerr=0;

    for(unsigned i=0; i<(1u<<31); i++){
        ap_fixed_vec3<1+16,1> x;
        
        ap_fixed<3+32,3> sum=0;
        for(int d=0; d<3; d++){
            x[d].setBits( urng() & 0x1FFFFul );
            sum += x[d] * x[d];
        }

        if(sum < 1){
            ap_fixed_vec3<1+16,1> xn;
            ap_ufixed<16,0> r;
            normalise_vector_16_v2(
                xn, r,
                x, sum
            );

            double true_r=sqrt(sum.to_double());
            double got_r=r.to_double();
            double err_r=std::abs(true_r-got_r);
            if(err_r > max_err){
                for(int d=0; d<3; d++){
                    std::cout<<" "<<x[d];
                }
                std::cout<<"   r^2="<<sum<<", got_r="<<r<<", true_r="<<true_r<<"\n";

                std::cout<<", err_r="<<err_r<<", log2="<<log2(err_r)<<"\n";
                max_err=err_r;
            }

            double cerr_r = 0;
            for(int d=0; d<3; d++){
                double true_x = x[d].to_double() / true_r;
                double got_x = xn[d].to_double();
                cerr_r = std::max(cerr_r, std::abs(true_x-got_x));
            }
            if(cerr_r > max_cerr){
                std::cout<<" cerr="<<cerr_r<<"\n";
                for(int d=0; d<3; d++){
                    double true_x = x[d].to_double() / true_r;
                    double got_x = xn[d].to_double();

                    std::cout<<"  "<<x[d]<<" -> " <<true_x<<" / "<<got_x;
                }
                std::cout<<"\n";
                max_cerr=cerr_r;
            }
        }
    }


}