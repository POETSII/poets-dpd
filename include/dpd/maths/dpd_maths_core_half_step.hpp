#ifndef dpd_maths_core_half_step_hpp
#define dpd_maths_core_half_step_hpp

#include "dpd/maths/dpd_maths_core.hpp"

#include "dpd/core/hash.hpp"

#include "dpd/core/logging.hpp"

#ifndef PDPD_TINSEL
#include <iostream>
#endif

namespace dpd_maths_core_half_step
{

    constexpr static double kT=  dpd_maths_core::kT;
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
            assert(0<=x[i] && x[i] < dims[i]);
        }
#ifndef PDPD_TINSEL
        //std::cerr<<"  alt: x="<<x<<", b.x'="<<b.x<<", v="<<b.v<<"\n";
#endif
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
    bool EnableLogging,
    class TScalar, class TVector,
    class TConservativeMap, class TDissipativeMap,
    class TBead1, class TBead2, class TForce
>
void calc_force(
    TScalar scale_inv_sqrt_dt,  // sqrt(24 * kT / dt)
    const TConservativeMap &conservative,
    const TDissipativeMap &sqrt_dissipative,

    uint64_t t_hash,

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
        
    auto conStrength=(TScalar)conservative(home_bead_type, other_bead_type);
    auto sqrtDissStrength=(TScalar)sqrt_dissipative(home_bead_type, other_bead_type); // This might be a constant

    auto dv = home.v - other.v;

    //std::cerr<<"  dv="<<dv<<"\n";

    TScalar wr = (TScalar(1) - dr);
        
    TScalar conForce = conStrength*wr;

    TVector dx_orig=dx;

    dx = dx * inv_dr;
        
    TScalar rdotv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    TScalar sqrt_gammap = sqrtDissStrength*wr;

    TScalar dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
    TScalar u = dpd_maths_core::default_hash(t_hash, home.get_hash_code(), other.get_hash_code());
    TScalar randScale = sqrt_gammap * scale_inv_sqrt_dt;
    TScalar randForce = randScale * u;

    TScalar dr0=r0-dr;
    TScalar hookeanForce=kappa*dr0;

#ifndef PDPD_TINSEL
    //std::cerr<<"  fhook="<<dx*hookeanForce<<", dx="<<dx<<",  xh="<<home.x<<", xt="<<other.x<<"\n";
    #endif

    TScalar scaled_force = conForce + dissForce + randForce + hookeanForce;

#ifndef PDPD_TINSEL
    //std::cerr<<"  "<<home.get_hash_code()<<" -> "<<other.get_hash_code()<<" : "<<conForce<<", "<<dissForce<<", "<<randForce<<", "<<hookeanForce<<"\n";
#endif

    force_home = dx * scaled_force;

    if(EnableLogging){
        if(ForceLogging::logger()){            
            double ddx[3]={dx_orig[0],dx_orig[1],dx_orig[2]};
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dx", 3,ddx);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dr", 1,&dr);
            double tt=scale_inv_sqrt_dt;
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-invrootdt", 1, &tt);
            double gammap=sqrt_gammap*sqrt_gammap;
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-gammap", 1, &gammap);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-rng", 1, &u);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-con",1, &conForce);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-diss", 1,&dissForce);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-rng-scale",1, &randScale);
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"dpd-rand",1, &randForce);
            double dpd_force=conForce + dissForce + randForce;
            double ff[3]={(double)dx[0]*dpd_force,(double)dx[1]*dpd_force,(double)dx[2]*dpd_force};
            ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"f_next_dpd", 3,ff);
            if(kappa!=0){
                double fh[3]={dx[0]*hookeanForce,dx[1]*hookeanForce,dx[2]*hookeanForce};
                ForceLogging::logger()->LogBeadPairProperty(home.get_hash_code(),other.get_hash_code(),"f_next_hookean", 3,fh);
            }
        }
    }

    //std::cerr<<"ref :   dr="<<dr<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<"\n";
        
}

template<
    bool EnableLogging,
    class TScalar, class TVector,
    class TForce
>
void calc_hookean_force(
    TScalar kappa, TScalar r0,
    TVector dx, TScalar dr,

    TForce & head_force, TForce &tail_force
) {
    TScalar dr0=r0-dr;
    TScalar hookeanForce=kappa*dr0; 

    head_force = - dx * (hookeanForce / dr);
    tail_force = - head_force;
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

    //std::cerr<<"  dx01="<<dx01<<", dx12="<<dx12<<", r0="<<r01<<", r1="<<r12<<"\n";

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
        if(absolute(b1Dotb2) > TScalar(0.000001) && sin_theta0>0){
            TScalar InvPrefactor = recip_pow_half(recip(cosPhiSq) - 1.0);

            // Add the restoring force if there is a preferred angle
            forceMag = forceMag*(cos_theta0 - sin_theta0*InvPrefactor);
        }

        //std::cerr<<"  kappa="<<kappa<<", cosPhiSq="<<cosPhiSq<<", forceMag="<<forceMag<<", maxProduct="<<magProduct<<"\n";

        //std::cerr<<"Dur: forceMag="<<forceMag<<"\n";

        headForce =((dx01*b1b2Overb1Sq)-dx12) * forceMag;
        assert(isfinite(headForce));

        tailForce = (dx01- (dx12*b1b2Overb2Sq)) * forceMag;
        assert(isfinite(tailForce));

        middleForce = -( headForce + tailForce );
    }else{
        headForce.clear();
        tailForce.clear();
        middleForce.clear();
    }
}

};

#endif