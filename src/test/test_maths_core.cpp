#include "dpd/core/vec3.hpp"
#include "dpd/core/dpd_state.hpp"
#include "dpd/maths/dpd_maths_core.hpp"
#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include <random>
#include <iostream>

struct input
{
    double dt;
    double con_strength, diss_strength;
    uint64_t t_hash;

    double kappa;
    double r0;


    struct{
        uint32_t hash;
        vec3r_t x;
        vec3r_t v;
        vec3r_t f;

        uint32_t get_bead_type() const
        { return bead_hash_get_bead_type(hash); }

        uint32_t get_hash_code() const
        { return hash; }
    }beads[2];
};

input random_input(std::mt19937_64 &rng)
{
    std::uniform_real_distribution<> udist;
    auto u=[&]() -> double { return udist(rng); };

    input i;
    i.dt=u()*0.2 + 0.01;
    i.con_strength = u();
    i.diss_strength = u();
    i.t_hash=rng();
    i.kappa=std::max(0.0, (u()-0.5)*10);
    i.r0=u()*0.4+0.3;
    i.beads[0].hash=rng()&0xFFFFul;
    i.beads[1].hash=rng()&0xFFFFul;
    for(int j=0; j<3; j++){
        i.beads[0].x[j]=u();
        i.beads[1].x[j]=u();
        i.beads[0].v[j]=u()*2-1;
        i.beads[1].v[j]=u()*2-1;
        i.beads[0].f[j]=u()*2-1;
        i.beads[1].f[j]=u()*2-1;
    }
    return i;
}

void test_input(const input &i)
{
    vec3r_t f_std, f_half;

    vec3r_t dx=i.beads[0].x-i.beads[1].x;
    double dr=dx.l2_norm();

    if( MIN_DISTANCE_CUTOFF_SQR < dr && dr <= 1 ){
        double scaled_inv_root_dt=pow_half(24*dpd_maths_core::kT / i.dt);
        dpd_maths_core::calc_force<double,vec3r_t>(
            0.5*i.dt, scaled_inv_root_dt,
            [&](unsigned, unsigned){ return i.con_strength; },
            [&](unsigned, unsigned){ return i.diss_strength; },
            i.t_hash,
            dx, dr,
            i.beads[0],
            i.beads[1],
            f_std
        );
        vec3r_t ff;
        dpd_maths_core::calc_hookean_force(
            i.kappa, i.r0, dx, dr, 1.0/dr, ff
        );
        f_std = f_std-ff;

        input ih=i;
        dpd_maths_core_half_step::update_mom(ih.dt, ih.beads[0]);
        dpd_maths_core_half_step::update_mom(ih.dt, ih.beads[1]);

        dpd_maths_core_half_step::calc_force<double,vec3r_t>(
            scaled_inv_root_dt,
            [&](unsigned, unsigned){ return i.con_strength; },
            [&](unsigned, unsigned){ return sqrt(i.diss_strength); },
            i.t_hash,
            dx, dr,
            i.kappa, i.r0,
            ih.beads[0],
            ih.beads[1],
            f_half
        );

        double err=(f_std-f_half).l2_norm();
        std::cerr<<"err="<<err<<"\n";
    }
}

int main()
{
    std::mt19937_64 rng;

    for(unsigned i=0; i<100;i++){
        auto in=random_input(rng);
        test_input(in);
    }
}
