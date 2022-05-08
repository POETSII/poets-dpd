#include "dpd/core/dpd_state_builder.hpp"
#include "dpd/core/dpd_state_io.hpp"

#include <random>
#include <regex>
#include <cmath>

int main(int argc, const char *argv[])
{
    double w=4, h=4, d=4;
    double density=3;
    double dt=0.05;

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

    const double g = 1.22074408460575947536;
    const double a1 = 1.0/g;
    const double a2 = 1.0/(g*g);
    const double a3 = 1.0/(g*g*g);

    builder.add_bead_type("W");
    builder.set_interaction_strength("W", "W", 20.0, 1.0);
    builder.add_polymer_type("P", {"W"}, {}, {});

    auto round_vec=[](const vec3r_t &x)
    {
        return vec3r_t{
            ldexp(round(ldexp(x[0],10)),-10),
            ldexp(round(ldexp(x[1],10)),-10),
            ldexp(round(ldexp(x[2],10)),-10)
        };
    };

    for(unsigned i=1; i<=n; i++){
        // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
        vec3r_t x{
            fmod(0.5+i*a1, 1.0)*w,
            fmod(0.5+i*a2, 1.0)*h,
            fmod(0.5+i*a3, 1.0)*d
        };
        vec3r_t x0=vec_wrap(round_vec(x), dims);
        builder.add_polymer("P", {{x}}, true);
    }

    state=builder.extract();
    state.dt=dt;

    write_world_state(std::cout, state);
}
