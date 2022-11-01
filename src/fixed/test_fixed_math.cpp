#include "fixed_math.hpp"
#include "fixed_math_v2.hpp"

#include "sliding_window.hpp"

#include <random>
#include <iostream>

/*
int main()
{
    std::mt19937_64 urng;

    double max_err=0;

    for(int i=0; i<1000000; i++){
        ap_fixed_vec3<1+16,1> x;
        
        ap_fixed<3+32,3> sum=0;
        for(int d=0; d<3; d++){
            x[d].setBits( urng() & 0x1FFFFul );
            sum += x[d] * x[d];
        }

        if(sum < 1){
            ap_fixed_vec3<1+16,1> xn;
            ap_ufixed<16,0> r;
            normalise_vector_16(
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
        }
    }


}
*/
