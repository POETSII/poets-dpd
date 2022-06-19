#include "dpd/core/dpd_state_builder.hpp"
#include "dpd/core/dpd_state_io.hpp"

#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb.hpp"

#include <random>
#include <regex>
#include <cmath>

int main(int argc, const char *argv[])
{
    double w=4, h=4, d=4;
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
    /*
    block([d],
        d : 6,
        allroots(x^(d+1) = x+1)
    );
    [x=0.9008649519489099*%i+0.6170934777839656,x=0.6170934777839656-0.9008649519489099*%i,x=0.2628696458512362*%i-0.8098578005941353,x=-0.2628696458512362*%i-0.8098578005941353,x=0.9525611952610316*%i-0.3636235193291831,x=-0.9525611952610316*%i-0.3636235193291831,x=1.112775684278706]
    */
    const double g = 1.112775684278706;
    const double a1 = 1.0/g;
    const double a2 = 1.0/(g*g);
    const double a3 = 1.0/(g*g*g);
    const double a4 = 1.0/(g*g*g*g);
    const double a5 = 1.0/(g*g*g*g*g);
    const double a6 = 1.0/(g*g*g*g*g*g);

    double r0=0.4;

    builder.add_bead_type("A");
    builder.add_bead_type("B");
    builder.add_bead_type("C");
    builder.set_interaction_strength("A", "A", 20.0, 1.0);
    builder.set_interaction_strength("A", "B", 20.0, 1.0);
    builder.set_interaction_strength("A", "C", 30.0, 1.0);

    builder.set_interaction_strength("B", "B", 20.0, 1.0);
    builder.set_interaction_strength("B", "C", 20.0, 1.0);

    builder.set_interaction_strength("C", "C", 20.0, 1.0);

    builder.add_polymer_type("P", {"A","B","C"}, {{40, r0, 0,1},{40, r0, 1,2}}, {});

    auto round_vec=[](const vec3r_t &x)
    {
        return vec3r_t{
            ldexp(round(ldexp(x[0],10)),-10),
            ldexp(round(ldexp(x[1],10)),-10),
            ldexp(round(ldexp(x[2],10)),-10)
        };
    };

    /* If a bead is going at 1/dt speed then it can
        cross an entire cell on one step, which breaks
        all kinds of things and beads will be lost. So
        we set a maximum speed of a quarter that.
    */
    double max_v = 0.25/dt;

    // Create an overall drift speed in some strange direction
    // that causes lots of wrapping and migrating.
    double drift_v = max_v / 3;
    vec3r_t drift={0.1,0.2,0.3};
    drift = drift * (drift_v/drift.l2_norm());

    double bead_v = max_v - drift_v;

    for(unsigned i=1; i<=n/3; i++){
        
        vec3r_t x{
            fmod(0.5+i*a1, 1.0)*w,
            fmod(0.5+i*a2, 1.0)*h,
            fmod(0.5+i*a3, 1.0)*d
        };

        vec3r_t dx{
            fmod(0.5+i*a4, 1.0)*2 - 1,
            fmod(0.5+i*a5, 1.0)*2 - 1,
            fmod(0.5+i*a6, 1.0)*2 - 1,
        };
        vec3r_t x0=vec_wrap(round_vec(x+dx*r0), dims);
        vec3r_t x1=vec_wrap(round_vec(x), dims);
        vec3r_t x2=vec_wrap(round_vec(x-dx*r0), dims);

        builder.add_polymer("P", {{x0},{x1},{x2}});
    }

    state=builder.extract();
    
    fprintf(stderr, "create_state_dimers : relaxing 100 steps at dt/10 using TBB.\n");
    state.dt=dt/10;
    NaiveDPDEngineHalfStepTBB engine;

    engine.Attach(&state);
    engine.Run(100);
    engine.Attach(nullptr);

    fprintf(stderr, "create_state_dimers : writing state.\n");

    state.dt=dt;
    write_world_state(std::cout, state);
}
