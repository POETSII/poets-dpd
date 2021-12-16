#include "dpd/engines/basic/basic_dpd_engine_v5_raw_orch.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/struct_to_c.hpp"

#include <random>
#include <regex>
#include <cmath>

int main(int argc, const char *argv[])
{
    WorldState state;

    vec3r_t dims{8,8,8};

    WorldStateBuilder builder(dims);

    double density=3;
    double w=dims[0], h=dims[1], d=dims[2];
    double volume=w*h*d;
    unsigned n=volume*density;

    const double g = 1.22074408460575947536;
    const double a1 = 1.0/g;
    const double a2 = 1.0/(g*g);
    const double a3 = 1.0/(g*g*g);

    builder.add_bead_type("A");
    builder.add_bead_type("B");
    builder.set_interaction_strength("A", "A", 20.0, 1.0);
    builder.set_interaction_strength("B", "B", 20.0, 1.0);
    builder.set_interaction_strength("A", "B", 30.0, 1.0);
    builder.add_polymer_type("P", {"A","B"}, {{10,1,0,1}}, {});

    auto round_vec=[](const vec3r_t &x)
    {
        return vec3r_t{
            ldexp(round(ldexp(x[0],10)),-10),
            ldexp(round(ldexp(x[1],10)),-10),
            ldexp(round(ldexp(x[2],10)),-10)
        };
    };

    for(unsigned i=1; i<=n/2; i++){
        // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
        vec3r_t x{
            fmod(0.5+i*a1, 1.0)*w,
            fmod(0.5+i*a2, 1.0)*h,
            fmod(0.5+i*a3, 1.0)*d
        };
        vec3r_t dx{
            fmod(0.1+i*a1, 1.0),
            fmod(0.1+i*a2, 1.0),
            fmod(0.1+i*a3, 1.0)
        };
        vec3r_t x0=vec_wrap(round_vec(x+dx*0.3), dims);
        vec3r_t x1=vec_wrap(round_vec(x-dx*0.3), dims);
        builder.add_polymer("P", {{x0},{x1}}, true);
    }

    state=builder.extract();

    BasicDPDEngineV5RawOrch engine;
    engine.Attach(&state);

    unsigned numSteps=100;
    engine.write_xml(std::cout, numSteps);
}
