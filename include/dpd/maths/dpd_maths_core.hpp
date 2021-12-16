#ifndef dpd_maths_core_hpp
#define dpd_maths_core_hpp

#include <cstdint>
#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"

#include "dpd/core/dpd_state.hpp"

#include "dpd/core/dpd_maths_primitives.hpp"

#ifndef TINSEL
#include <iostream>
#endif

namespace dpd_maths_core
{

    static constexpr double kT = 1.0;

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
    inline int32_t uint32_to_int32(uint32_t x)
    {
        // Sigh, Avoid undefined behaviour. Compiler should optimise it out.
        // https://stackoverflow.com/a/13208789
        int32_t res;
        if (x <= INT32_MAX) {
            res=static_cast<uint32_t>(x);
        }else{
            assert(x >= (uint32_t)INT32_MIN);
            res= static_cast<uint32_t>(x - INT32_MIN) + INT32_MIN;
        }
        assert(x == uint32_t(res) );  // int32_t -> uint32_t is well-defined
        return res;
    }
};

inline float default_hash(uint64_t base, uint32_t s1, uint32_t s2)
{
    uint32_t ru=hash_rng_sym(base, s1, s2);
    int32_t rs=detail::uint32_to_int32(ru);  // in [-2^31,2^31)
    //const double scale=ldexp(1.0,-32) / sqrt(1/3.0); // gives stddev of 1 (same as groot-warren paper)
    const float scale=0.00000000023283064365386962890625f; // Gives range of [-0.5,0.5]  (same as Osprey-DPD)
    float u = rs * scale; 
    /*
    static double u_sum_sqr=0, u_sum=0;
    static unsigned u_count=0;
    u_sum_sqr += u*u;
    u_sum += u;
    u_count += 1;
    if((u_count % 10000)==0){
        std::cerr<<"ucount="<<u_count<<", mean="<<u_sum/u_count<<", std="<<sqrt(u_sum_sqr/u_count)<<"\n";
    }
    */
    return u;
}

inline float default_hash(uint64_t base, const BeadHash &s1, const BeadHash &s2)
{
    return default_hash(base, s1.hash, s2.hash);
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
    TForce &f,
    TBead &b
) {
    // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

    b.v += (b.f + f) * half(dt);
    b.f = f;
    f.clear();
}

template<
    class TScalar, class TVector,
    class TConservativeMap, class TDissipativeMap,
    class TBead1, class TBead2, class TForce
>
void calc_force(
    TScalar lambda_dt,
    TScalar scaled_inv_sqrt_dt, // sqrt(24 * kT / dt)
    const TConservativeMap &conservative,
    const TDissipativeMap &dissipative,

    uint64_t t_hash,

    TVector dx, TScalar dr,

    const TBead1 &home,
    const TBead2 &other,

    TForce & force_home
) {
    assert(&home != &other);
    assert(home.get_hash_code() != other.get_hash_code());
    assert(dr < 1);

    TScalar inv_dr = recip(dr);

    auto home_bead_type=home.get_bead_type();
    auto other_bead_type=other.get_bead_type();
        
    auto conStrength=conservative(home_bead_type, other_bead_type);
    auto dissStrength=dissipative(home_bead_type, other_bead_type); // This might be a constant

    vec3r_t hb_v = home.v + home.f * lambda_dt;
    vec3r_t ob_v = other.v + other.f * lambda_dt;
    vec3r_t dv = hb_v - ob_v;

    //std::cerr<<"  dv="<<dv<<"\n";

    TScalar wr = (1 - dr);
    TScalar wr2 = wr*wr;
        
    TScalar conForce = conStrength*wr;
        
    TScalar rdotv = (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
    TScalar gammap = dissStrength*wr2;

    TScalar dissForce = -gammap*rdotv;
    TScalar u = default_hash(t_hash, home.get_hash_code(), other.get_hash_code());
    TScalar randForce = pow_half(gammap) * scaled_inv_sqrt_dt * u;

    TScalar scaled_force = (conForce + dissForce + randForce) * inv_dr;

    force_home = dx * scaled_force;
}

template<class TScalar, class TVector, class TForce>
void calc_hookean_force(
    TScalar kappa,
    TScalar r0,
    const TVector &dx, // Distance from head to tail
    TScalar r,         // == dx.l2_norm()
    TScalar inv_r, // Quite possibly already calculated
    TForce &force_head
){
    // The force scale is just kappa*(r-r0).
    TScalar dr=(r-r0);
    TScalar force=kappa*dr;

    // Division by r is just to get (dx/r)
    force_head=dx * (force * inv_r);
    //std::cerr<<"  dx="<<dx<<", dr="<<r<<", r0="<<r0<<"\n";
    //std::cerr<<"  fhead="<<force_head<<"\n";
}

template<class TScalar, class TVector, class TForce>
void calc_angle_force(
    TScalar kappa,
    TScalar cos_theta0,
    TScalar sin_theta0,
    const TVector &dx01, TScalar r01,
    const TVector &dx12, TScalar r12,
    TForce &headForce,
    TForce &middleForce,
    TForce &tailForce
){
    TScalar magProduct = r01*r12;

    if(magProduct > TScalar(0.0001))
    {
        TScalar b1MagSq		= r01*r01;
        TScalar b2MagSq		= r12*r12;
        TScalar b1Dotb2		= dx01[0]*dx12[0] + dx01[1]*dx12[1] + dx01[2]*dx12[2];
        TScalar b1b2Overb1Sq	= b1Dotb2/b1MagSq;
        TScalar b1b2Overb2Sq	= b1Dotb2/b2MagSq;
        TScalar cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

        TScalar forceMag = kappa/magProduct;

        // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0
        if(absolute(b1Dotb2) > TScalar(0.000001) && sin_theta0>0){
            TScalar InvPrefactor = recip_pow_half(recip(cosPhiSq) - 1);

            // Add the restoring force if there is a preferred angle
            forceMag = forceMag*(cos_theta0 - sin_theta0*InvPrefactor);
        }

        headForce =((dx01*b1b2Overb1Sq)-dx12) * forceMag;
        assert(isfinite(headForce));

        tailForce = (dx01- (dx12*b1b2Overb2Sq)) * forceMag;
        assert(isfinite(tailForce));

        middleForce = -( headForce + tailForce );
    }
}

};

#endif