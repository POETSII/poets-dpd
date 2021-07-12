#ifndef dpd_maths_core_hpp
#define dpd_maths_core_hpp

#include <cstdint>
#include "vec3.hpp"

namespace dpd_maths_core
{

    double half(double x)
    { return x*0.5; }

    float half(float x)
    { return x*0.5f; }

struct bead_concept_t
{
    using vec_t = vec3r_t;

    uint32_t get_bead_id() const;
    uint32_t get_polymer_id() const;

    unsigned get_bead_type() const;
    unsigned get_polymer_type() const;
};

namespace detail
{
    inline uint64_t splitmix64(uint64_t z)
    {
        z = (z ^ (z >> 30)) * uint64_t(0xBF58476D1CE4E5B9ull);
        z = (z ^ (z >> 27)) * uint64_t(0x94D049BB133111EBull);
        return z ^ (z >> 31);
    }

    inline uint64_t RandU64(uint64_t base, uint32_t s1, uint32_t s2) 
    {
        uint64_t z = base + ( ( uint64_t(std::max(s1,s2))<<32) | std::min(s1,s2) );
        return splitmix64(z);
    }
};

inline double default_hash(uint64_t base, uint32_t s1, uint32_t s2)
{
    uint32_t ru=(uint32_t)detail::RandU64(base,s1,s2);
    const double scale=ldexp(2.0, -32);
    double u = ru * scale; 
    return u-0.5;  // Gives range of [-0.5,0.5]  (same as Osprey-DPD)
}

template<class TScalar, class TDims, class TBead>
void update_pos(
    TScalar dt,
    const TDims &dims,
    TBead &b
)
{
    auto x = b.x + b.v*dt + b.f*half(dt*dt);
    for(int i=0; i<3; i++){
        x[i] += (x[i]<0 ? dims[i] : 0);
        x[i] -= (x[i]>=dims[i] ? dims[i] : 0);
    }
    b.x = x;
}

template<class TScalar, class TBead, class TForce>
void update_mom(
    TScalar dt,
    const TForce &f,
    TBead &b
) {
    // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

    b.v += (b.f + f) * (dt * 0.5);
    b.f = f;
}


//! Returns true if there is an interaction, with force_home the force on bead home.
template<
    class TScalar,
    class TConservativeMap, class TDissipativeMap, class TRandHash,
    class TBead1, class TBead2, class TDistance, class TForce
>
bool calc_force(
    TScalar lambda_dt,
    TScalar inv_sqrt_dt,
    const TConservativeMap &conservative,
    const TDissipativeMap &dissipative,

    uint64_t hash_base,
    const TRandHash &hash_unif,  // maps (uint64_t,uint32_t,uint32_t) -> TScalar. Uniform in [-0.5,0.5]

    const TBead1 &home,
    const TBead2 &other,

    const TDistance &wrap_delta,  // Vector offset to handle wrapping

    TForce & force_home
) {
    assert(&home != &other);
    assert(home.get_bead_id() != other.get_bead_id());

    vec3r_t dx = vec3r_t(home.x) - vec3r_t(other.x) - wrap_delta; 

    TScalar dr=dx.l2_norm();
    TScalar inv_dr=1.0/dr;

    if(dr>=1 || dr < 0.000000001){
        return false;
    }

    auto home_bead_type=home.get_bead_type();
    auto other_bead_type=other.get_bead_type();
        
    auto conStrength=conservative(home_bead_type, other_bead_type);
    auto dissStrength=dissipative(home_bead_type, other_bead_type); // This might be a constant

    vec3r_t hb_v = home.v + home.f * lambda_dt;
    vec3r_t ob_v = other.v + other.f * lambda_dt;
    vec3r_t dv = hb_v - ob_v;

    TScalar wr = (1.0 - dr);
    TScalar wr2 = wr*wr;
        
    TScalar conForce = conStrength*wr;
        
    TScalar rdotv = (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
    TScalar gammap = dissStrength*wr2;

    TScalar dissForce = -gammap*rdotv;
    TScalar u = hash_unif(hash_base, home.get_bead_id(), other.get_bead_id());
    TScalar randForce = sqrt(gammap) * inv_sqrt_dt * u;

    TScalar scaled_force = (conForce + dissForce + randForce) * inv_dr;

    force_home = dx * scaled_force;

    return true;
}

};

#endif