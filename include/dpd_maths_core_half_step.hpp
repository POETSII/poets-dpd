#ifndef dpd_maths_core_half_step_hpp
#define dpd_maths_core_half_step_hpp

#include "dpd_maths_core.hpp"

#include "hash.hpp"

namespace dpd_maths_core_half_step
{

    using dpd_maths_core::half;
    using dpd_maths_core::recip;
    using dpd_maths_core::recip_sqrt;

    using dpd_maths_core::default_hash;

    // After this method b.f is dead
    template<class TScalar, class TDims, class TBead>
    void update_pos(
        TScalar dt,
        const TDims &dims,
        TBead &b
    )
    {
        // x(t+dt) = x(t) + v(t)*dt + f(t)*dt*dt
        // v(t+dt/2) = v(t) + dt*f(t)/2

        auto x = b.x + b.v*dt + b.f*half(dt*dt);
        for(int i=0; i<3; i++){
            x[i] += (x[i]<0 ? dims[i] : 0);
            x[i] -= (x[i]>=dims[i] ? dims[i] : 0);
        }
        b.x = x;
        b.v = b.v + b.f * half(dt);
        b.f.clear();
    }

template<class TScalar, class TBead>
void update_mom(
    TScalar dt,
    TBead &b
) {
    // v(t+dt) = v(t+dt/2) + dt*f(t+dt)/2
    // v(t+dt) = v(t) + dt*f(t)/2 + dt*f(t+dt)/2
    // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

    b.v += b.f * half(dt);
}

template<
    class TScalar, class TVector,
    class TConservativeMap, class TDissipativeMap,
    class TBead1, class TBead2, class TForce
>
void calc_force(
    TScalar inv_sqrt_dt,
    const TConservativeMap &conservative,
    const TDissipativeMap &dissipative,

    uint32_t t_hash,

    TVector dx, TScalar dr,
    TScalar kappa, TScalar r0,

    const TBead1 &home,
    const TBead2 &other,

    TForce & force_home
) {
    assert(&home != &other);
    assert(home.get_hash_code() != other.get_hash_code());
    assert(0 < dr && dr < 1);

    TScalar inv_dr = recip(dr);

    auto home_bead_type=home.get_bead_type();
    auto other_bead_type=other.get_bead_type();
        
    auto conStrength=conservative(home_bead_type, other_bead_type);
    auto dissStrength=dissipative(home_bead_type, other_bead_type); // This might be a constant

    auto dv = home.v - other.v;

    TScalar wr = (1.0 - dr);
    TScalar wr2 = wr*wr;
        
    TScalar conForce = conStrength*wr;
        
    TScalar rdotv = (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
    TScalar gammap = dissStrength*wr2;

    TScalar dissForce = -gammap*rdotv;
    TScalar u = dpd_maths_core::default_hash(t_hash, home.get_hash_code(), other.get_hash_code());
    TScalar randForce = sqrt(gammap) * inv_sqrt_dt * u;

    TScalar dr0=r0-dr;
    TScalar hookeanForce=kappa*dr0;

    TScalar scaled_force = (conForce + dissForce + randForce + hookeanForce) * inv_dr;

    force_home = dx * scaled_force;

    //std::cerr<<"ref :   dr="<<dr<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<"\n";
        
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

    if(magProduct > 0.0001)
    {
        TScalar b1MagSq		= r01*r01;
        TScalar b2MagSq		= r12*r12;
        TScalar b1Dotb2		= dx01[0]*dx12[0] + dx01[1]*dx12[1] + dx01[2]*dx12[2];
        TScalar b1b2Overb1Sq	= b1Dotb2/b1MagSq;
        TScalar b1b2Overb2Sq	= b1Dotb2/b2MagSq;
        TScalar cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

        TScalar forceMag = kappa/magProduct;

        // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0
        if(fabs(b1Dotb2) > 0.000001 && sin_theta0>0){
            TScalar InvPrefactor = recip_sqrt(recip(cosPhiSq) - 1.0);

            // Add the restoring force if there is a preferred angle
            forceMag = forceMag*(cos_theta0 - sin_theta0*InvPrefactor);
        }

        //std::cerr<<"Dur: forceMag="<<forceMag<<"\n";

        headForce =((dx01*b1b2Overb1Sq)-dx12) * forceMag;
        assert(isfinite(headForce));

        tailForce = (dx01- (dx12*b1b2Overb2Sq)) * forceMag;
        assert(isfinite(tailForce));

        middleForce = -( headForce + tailForce );
    }
}

};

#endif