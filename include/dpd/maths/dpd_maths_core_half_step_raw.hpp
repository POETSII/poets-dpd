#ifndef dpd_maths_core_half_step_raw_hpp
#define dpd_maths_core_half_step_raw_hpp

#include "dpd/maths/dpd_maths_core.hpp"

#include "dpd/core/hash.hpp"

#ifndef TINSEL
#include <iostream>
#endif

namespace dpd_maths_core_half_step_raw
{

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


        // auto x = b.x + b.v*dt + b.f*half(dt*dt);
        // auto x = b.x + dt * ( b.v + b.f*half(dt) );
        float x[3];
        vec3_mul(x, b.f, half(dt));
        vec3_add(x, b.v);
        vec3_mul(x, dt);
        vec3_add(x, b.x);

        
        for(int i=0; i<3; i++){
            x[i] += (x[i]<0 ? dims[i] : 0);
            x[i] -= (x[i]>=dims[i] ? dims[i] : 0);
            if(x[i] < 0){
                printf("<\n");
            }
            if(x[i] >= dims[i]){
                printf(">\n");
            }
        }
        vec3_copy(b.x, x);
        vec3_add_mul(b.v, b.f, half(dt));
        vec3_clear(b.f);
    }

    // Assumes that x is an integer with 22 fractional bits
    template<class TScalar, class TBead>
    void update_pos_f22(
        TScalar dt,
        const int32_t dims_f22[3],
        TBead &b
    )
    {
        // x(t+dt) = x(t) + v(t)*dt + f(t)*dt*dt
        // v(t+dt/2) = v(t) + dt*f(t)/2


        // auto x = b.x + b.v*dt + b.f*half(dt*dt);
        // auto x = b.x + dt * ( b.v + b.f*half(dt) );
        float dx[3];
        vec3_mul(dx, b.f, half(dt));
        vec3_add(dx, b.v);
        vec3_mul(dx, dt);

        int32_t orig[3];
        vec3_copy(orig, b.x_f22);

        int32_t nx[3];
        for(int i=0; i<3; i++){
            nx[i] = b.x_f22[i] + round_impl( dx[i] * float(1<<22) );
        }

        
        for(int i=0; i<3; i++){
            nx[i] += (nx[i]<0 ? dims_f22[i] : 0);
            nx[i] -= (nx[i]>=dims_f22[i] ? dims_f22[i] : 0);
            if(nx[i] < 0){
                printf("<\n");
            }
            if(nx[i] >= dims_f22[i]){
                printf(">\n");
            }
        }
        vec3_copy(b.x_f22, nx);
        
        vec3_add_mul(b.v, b.f, half(dt));
        vec3_clear(b.f);
    }

template<class TScalar, class TBead>
void update_mom(
    TScalar dt,
    TBead &b
) {
    // v(t+dt) = v(t+dt/2) + dt*f(t+dt)/2
    // v(t+dt) = v(t) + dt*f(t)/2 + dt*f(t+dt)/2
    // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

    vec3_add_mul(b.v, b.f, half(dt));
}

template<
    bool EnableLogging, class TScalar, class TVector, class TForce
>
void calc_force(
    TScalar scale_inv_sqrt_dt,
    uint64_t t_hash,

    TVector dx, TScalar dr,
    TScalar kappa, TScalar r0,

    TScalar conStrength,
    TScalar sqrtDissStrength,
    BeadHash home_hash,
    BeadHash other_hash,
    const TVector &home_v,
    const TVector &other_v,

    TForce & force_home
) {
    assert(home_hash != other_hash);
    assert(0 < dr && dr < 1);

    TScalar inv_dr = recip(dr);
        
    TVector dv;
    vec3_sub(dv, home_v, other_v);

    TScalar wr = (1 - dr);
        
    TScalar conForce = conStrength*wr;

    float dx_orig[3]={dx[0],dx[1],dx[2]};

    vec3_mul(dx, inv_dr);
        
    TScalar rdotv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    TScalar sqrt_gammap = sqrtDissStrength*wr;

    TScalar dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
    TScalar u = dpd_maths_core::default_hash(t_hash, home_hash, other_hash);
    TScalar randScale = sqrt_gammap * scale_inv_sqrt_dt ;
    TScalar randForce = randScale * u;

    TScalar dr0=r0-dr;
    TScalar hookeanForce=kappa*dr0;

    TScalar scaled_force = conForce + dissForce + randForce + hookeanForce;

#ifndef TINSEL
    //std::cerr<<"  DUT: home="<<home_hash.hash<<", other="<<other_hash.hash<<", t_hash="<<t_hash<<", dx=("<<(dx[0]*dr)<<","<<(dx[1]*dr)<<","<<(dx[2]*dr)<<"), r="<<dr<<", u="<<u<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<", hook="<<hookeanForce<<"\n";
    //std::cerr<<"     sqrt_gammap="<<sqrt_gammap<<", rdotv="<<rdotv<<", pow_half(dissStrength)="<<sqrtDissStrength<<"\n";
#endif

    vec3_mul(force_home, dx , scaled_force);

    for(int i=0; i<3; i++){
        assert(abs(force_home[i])<1000);
    }

    if(EnableLogging){
        if(ForceLogging::logger()){
           
            double ddx[3]={dx_orig[0],dx_orig[1],dx_orig[2]};
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dx", 3,ddx);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dr", 1,&dr);
            double dd=sqrtDissStrength*sqrtDissStrength;
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-diss-strength", 1, &dd);
            double tt=scale_inv_sqrt_dt;
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-invrootdt", 1, &tt);
            double gammap=sqrt_gammap*sqrt_gammap;
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-gammap", 1, &gammap);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-rng", 1, &u);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-con",1, &conForce);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-diss", 1,&dissForce);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-rng-scale",1, &randScale);
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"dpd-rand",1, &randForce);
            double dpd_force=conForce + dissForce + randForce;
            double ff[3]={dx[0]*dpd_force,dx[1]*dpd_force,dx[2]*dpd_force};
            ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"f_next_dpd", 3,ff);
            if(kappa!=0){
                double fh[3]={dx[0]*hookeanForce,dx[1]*hookeanForce,dx[2]*hookeanForce};
                ForceLogging::logger()->LogBeadPairProperty(home_hash,other_hash,"f_next_hookean", 3,fh);
            }
        }
    }
}

template<bool AlwaysStraight, class TScalar, class TVector, class TForce>
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

    if(magProduct > TScalar(0.0001))
    {
        TScalar b1MagSq		= r01*r01;
        TScalar b2MagSq		= r12*r12;
        TScalar b1Dotb2		= dx01[0]*dx12[0] + dx01[1]*dx12[1] + dx01[2]*dx12[2];
        TScalar b1b2Overb1Sq	= b1Dotb2/b1MagSq;
        TScalar b1b2Overb2Sq	= b1Dotb2/b2MagSq;
        TScalar cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

        

        TScalar forceMag = kappa/magProduct;

        if(!AlwaysStraight){
            // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0
            if(absolute(b1Dotb2) > TScalar(0.000001) && sin_theta0>0){
                TScalar InvPrefactor = recip_pow_half(recip(cosPhiSq) - TScalar(1));

                // Add the restoring force if there is a preferred angle
                forceMag = forceMag*(cos_theta0 - sin_theta0*InvPrefactor);
            }
        }else{
            assert(sin_theta0==0);
        }

        //std::cerr<<"  kappa="<<kappa<<", cosPhiSq="<<cosPhiSq<<", forceMag="<<forceMag<<", maxProduct="<<magProduct<<"\n";

        //std::cerr<<"Dur: forceMag="<<forceMag<<"\n";

        //headForce =((dx01*b1b2Overb1Sq)-dx12) * forceMag;
        vec3_mul(headForce, dx01, b1b2Overb1Sq );
        vec3_sub(headForce, dx12);
        vec3_mul(headForce, forceMag);
        assert(isfinite(headForce));
        for(int i=0; i<3; i++){
            assert(abs(headForce[i])<100);
        }

        //tailForce = (dx01- (dx12*b1b2Overb2Sq)) * forceMag;
        vec3_mul(tailForce, dx12, b1b2Overb2Sq );
        vec3_sub(tailForce, dx01, tailForce);
        vec3_mul(tailForce, forceMag);
        assert(isfinite(tailForce));
        for(int i=0; i<3; i++){
            assert(abs(tailForce[i])<100);
        }

        //middleForce = -( headForce + tailForce );
        vec3_add(middleForce, headForce, tailForce);
        vec3_neg(middleForce);
    }else{
        vec3_clear(headForce);
        vec3_clear(tailForce);
        vec3_clear(middleForce);
    }
}

};

#endif