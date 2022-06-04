#include "dpd/core/dpd_state_builder.hpp"
#include "dpd/core/dpd_state_io.hpp"

#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb.hpp"

#include <random>
#include <regex>
#include <cmath>

int main(int argc, const char *argv[])
{
    double w=128, h=128, d=8;
    double density=3;
    double dt=0.01;

    if(argc>1){
        w=std::atoi(argv[1]);
    }
    if(argc>2){
        h=std::atoi(argv[2]);
    }
    if(argc>3){
        d=std::atoi(argv[3]);
    }
    if(argc>4){
        density=std::atof(argv[4]);
    }
    if(argc>5){
        dt=std::atof(argv[5]);
    }

    WorldState state;

    vec3r_t dims{w,h,d};

    WorldStateBuilder builder(dims);

    double volume=w*h*d;
    unsigned n=round(volume*density);

    // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const double g = 1.22074408460575947536;
    const double a1 = 1.0/g;
    const double a2 = 1.0/(g*g);
    const double a3 = 1.0/(g*g*g);


    builder.add_bead_type("W");
    builder.add_bead_type("A");
    builder.add_bead_type("B");

    builder.set_interaction_strength("W", "W", 50.0, 5.0);
    builder.set_interaction_strength("W", "A", 50.0, 5.0);
    builder.set_interaction_strength("W", "B", 50.0, 5.0);

    builder.set_interaction_strength("A", "B", "10+ux*20+uy*20", 5.0);

    builder.set_interaction_strength("A", "A", "10+ux*40", 5.0);
    builder.set_interaction_strength("B", "B", "10+uy*40", 5.0);

    builder.add_polymer_type("W", {"W"}, {}, {});
    builder.add_polymer_type("A", {"A"}, {}, {});
    builder.add_polymer_type("B", {"B"}, {}, {});


    auto round_vec=[](const vec3r_t &x)
    {
        return vec3r_t{
            ldexp(round(ldexp(x[0],10)),-10),
            ldexp(round(ldexp(x[1],10)),-10),
            ldexp(round(ldexp(x[2],10)),-10)
        };
    };

    std::mt19937_64 urng;
    std::uniform_real_distribution<> udist;

    for(unsigned i=1; i<=n; i++){
        
        vec3r_t x{
            fmod(0.5+i*a1, 1.0)*w,
            fmod(0.5+i*a2, 1.0)*h,
            fmod(0.5+i*a3, 1.0)*d
        };
        double u=udist(urng);
        if(u<0.2){
            builder.add_polymer("A", {{x}});
        }else if(u<0.4){
            builder.add_polymer("B", {{x}});
        }else{
            builder.add_polymer("W", {{x}});
        }
    }

    state=builder.extract();
    state.dt=dt;
    write_world_state(std::cout, state);
}
