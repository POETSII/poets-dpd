#include "dpd/core/cvec3.hpp"

#include <random>
#include <iostream>

void require(bool cond, const char *msg)
{
    if(!cond){
        throw std::runtime_error(msg);
    }
}

#define REQUIRE(cond) require(cond, #cond)

int main()
{
    {
        CVec3Half h;
    
        REQUIRE(h.get_x()==0);
        REQUIRE(h.get_y()==0);
        REQUIRE(h.get_z()==0);

        std::cerr<<"Testing integers\n";
        for(int x=-255; x<=255; x++){
            h=CVec3Half(x, 0, 0);
            //std::cerr<<"x="<<x<<", h.get_x()=="<<h.get_x()<<"\n";
            REQUIRE(h.get_x()==x);
            REQUIRE(h.get_y()==0);
            REQUIRE(h.get_z()==0);

            for(int y=-255; y<=255; y+=7){
                h=CVec3Half(x, y, 0);
                //std::cerr<<"x="<<x<<", h.get_x()=="<<h.get_x()<<", y="<<y<<", h.get_y()="<<h.get_y()<<"\n";
                REQUIRE(h.get_x()==x);
                REQUIRE(h.get_y()==y);
                REQUIRE(h.get_z()==0);
            
                for(int z=-255; z<=255; z+=13){
                    h=CVec3Half(x, y, z);
                    //std::cerr<<"x="<<x<<", h.get_x()=="<<h.get_x()<<", y="<<y<<", h.get_y()="<<h.get_y()<<"\n";
                    REQUIRE(h.get_x()==x);
                    REQUIRE(h.get_y()==y);
                    REQUIRE(h.get_z()==z);
                }
            }
        }

        std::mt19937_64 urng;
        std::uniform_real_distribution<> udist;

        std::cerr<<"Testing floats in [-100,100]^3\n";
        for(int i=0; i<1000000; i++){
            double x=(2*udist(urng)-1)*100;
            double y=(2*udist(urng)-1)*100;
            double z=(2*udist(urng)-1)*100;
            
            CVec3Half h(x,y,z);

            if(h.is_zero()){
                REQUIRE( std::abs(x) <= ldexp(0.5, -14));
                REQUIRE( std::abs(y) <= ldexp(0.5, -14));
                REQUIRE( std::abs(z) <= ldexp(0.5, -14));
            }else{
                double largest=std::max(std::abs(x),std::max(std::abs(y),std::abs(z)));

                //std::cerr<<"  x="<<x<<", get_x()="<<h.get_x()<<", largest="<<largest<<", err="<<std::abs((h.get_x() - x) / largest)<<"\n";
                REQUIRE( std::abs((h.get_x() - x) / largest) < 1.0/128 );
                REQUIRE( std::abs((h.get_y() - y) / largest) < 1.0/128 );
                REQUIRE( std::abs((h.get_z() - z) / largest) < 1.0/128 );

                double tr=vec3r_t(x,y,z).l2_norm();
                double gr=h.get_vec3r().l2_norm();
                //std::cerr<<" tr="<<tr<<", gr="<<gr<<"\n";
                REQUIRE( std::abs( (tr-gr) / gr ) < 1.0/64.0 ); // Not sure exactly what relative precision in magnitude should be.
            }
        }

        std::cerr<<"Testing floats closer to zero.\n";
        for(int i=0; i<1000000; i++){
            double x=(2*udist(urng)-1);
            double y=(2*udist(urng)-1);
            double z=(2*udist(urng)-1);

            x=x*x*x*100;
            y=y*y*y*100;
            z=z*z*z*100;
            
            CVec3Half h(x,y,z);

            if(h.is_zero()){
                REQUIRE( std::abs(x) <= ldexp(0.5, -14));
                REQUIRE( std::abs(y) <= ldexp(0.5, -14));
                REQUIRE( std::abs(z) <= ldexp(0.5, -14));
            }else{
                double largest=std::max(std::abs(x),std::max(std::abs(y),std::abs(z)));

                //std::cerr<<"  x="<<x<<", get_x()="<<h.get_x()<<", largest="<<largest<<", err="<<std::abs((h.get_x() - x) / largest)<<"\n";
                REQUIRE( std::abs((h.get_x() - x) / largest) < 1.0/128 );
                REQUIRE( std::abs((h.get_y() - y) / largest) < 1.0/128 );
                REQUIRE( std::abs((h.get_z() - z) / largest) < 1.0/128 );

                double tr=vec3r_t(x,y,z).l2_norm();
                double gr=h.get_vec3r().l2_norm();
                //std::cerr<<" tr="<<tr<<", gr="<<gr<<"\n";
                REQUIRE( std::abs( (tr-gr) / gr ) < 1.0/64.0 ); // Not sure exactly what relative precision in magnitude should be.
            }
        }
    }

    std::cerr<<"Ok\n";


}