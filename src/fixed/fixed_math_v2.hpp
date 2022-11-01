#ifndef fixed_math_v2_hpp
#define fixed_math_v2_hpp

#include "fixed_math.hpp"

struct fixed_maths_v2
{
    static const int BEAD_ID_BITS = 32;
    static const int ROUND_SEED_BITS = 32;
    static const int BEAD_TYPE_BITS = 4;

    static const int CONSERVATIVE_FRAC_BITS = 0; // This is always integer in practise, AFAICT
    static const int CONSERVATIVE_INT_BITS = 7;  // Never seen values more than 100
    using conservative_strength_t = ap_ufixed<CONSERVATIVE_INT_BITS+CONSERVATIVE_FRAC_BITS,CONSERVATIVE_FRAC_BITS>;

    // Dissipative strength is generally in the range [1,64) with fractional precision of 0.5
    // However, we want sqrt of that, which means it will turn into something irrational. So
    // assuming dissipative strength versus conservative is the important ratio, we need
    // reasonable precision here. So let's assume that sqrt([1,64))->[1,8) and we need decent
    // fractional precision
    static const int SQRT_DISSIPATIVE_FRAC_BITS = 16;
    static const int SQRT_DISSIPATIVE_INT_BITS = 3;
    using sqrt_dissipative_strength_t = ap_ufixed<SQRT_DISSIPATIVE_INT_BITS+SQRT_DISSIPATIVE_FRAC_BITS,SQRT_DISSIPATIVE_FRAC_BITS>;

    static const int POS_FRAC_BITS = 16;
    using pos_t = ap_fixed_vec3<POS_FRAC_BITS, 0>;
    
    // Position delta must be in (-1,+1)^3 for distance check to pass
    using pos_delta_t = ap_fixed_vec3<POS_FRAC_BITS+1, 1>;

    // Distance is always in (0,1) for distance check to pass
    using distance_sqr_t = ap_ufixed<2*POS_FRAC_BITS,0>;

    // Distance is always in (0,1) for distance check to pass
    using distance_t = ap_ufixed<POS_FRAC_BITS,0>;

    // This is the fractional accuracy of normalised pure directional vectors
    static const int DIRECTION_FRAC_BITS = 16;
    using direction_t = ap_fixed_vec3<1+DIRECTION_FRAC_BITS,1>;

    // Maximum allowable speed is 1/dt, or realistically 0.5/dt. For dt of 0.01
    // we get max dt of 0.5/0.01 = 50. I'm going to call that 31 for convenience.
    // Max relative passing speed is then 62.
    static const int VEL_INT_BITS = 5;
    static const int VEL_FRAC_BITS = 14;
    using velocity_part_t = ap_fixed<1+VEL_INT_BITS+VEL_FRAC_BITS,1+VEL_INT_BITS>;
    using velocity_t = ap_fixed_vec3<1+VEL_INT_BITS+VEL_FRAC_BITS,1+VEL_INT_BITS>;
    using velocity_delta_t = ap_fixed_vec3<1+1+VEL_INT_BITS+VEL_FRAC_BITS,1+VEL_INT_BITS>;

    // Fractional precision of the various force scalar components
    static const int COMPONENT_FRAC_BITS = 16;

    static const int FORCE_INT_BITS = 10;  // Should never be more than 1024 for relaxed system?
    static const int FORCE_FRAC_BITS = 16; // Seems reasonable...?
    using force_part_t = ap_fixed<1+FORCE_INT_BITS+FORCE_FRAC_BITS,FORCE_FRAC_BITS>;
    using force_t = ap_fixed_vec3<1+FORCE_INT_BITS+FORCE_FRAC_BITS,FORCE_FRAC_BITS>;

    struct bead_info_t
    {
        ap_uint<BEAD_ID_BITS> id;
        ap_uint<BEAD_TYPE_BITS> type;
        pos_t x;
        velocity_t v;
    };

    __attribute__((always_inline)) static void calc_force_fix(
        const bead_info_t &hb,
        const bead_info_t &ob,
        ap_fixed_vec3<2,1> ob_cell_delta,
        const ap_uint<ROUND_SEED_BITS> &round_seed,
        force_t &hb_f
    ){
        // TODO: Make run-time loadable
        static const conservative_strength_t con_table[1024]={1,2,3,1,2,3,1,2,3,0,2,1,3, 0.23211};
        static const sqrt_dissipative_strength_t sqrt_diss_table[1024]={1,2,3,1,1,2,3,1,2,3,0,2,1,3,0,2,0.1244,2,7,0.21232};

        ap_fixed_vec3<2+POS_FRAC_BITS,POS_FRAC_BITS> dx_full;
        ap_ufixed<3+POS_FRAC_BITS*2,POS_FRAC_BITS*2> dr2_full = 0;
        pos_delta_t dx;
#pragma HLS UNROLL
        for(int d=0; d<3; d++){
            auto dx_full = ob.x[d] - hb.x[d] - ob_cell_delta[d];
            dr2_full += dx_full * dx_full;
            dx[d] = dx_full;
        }
        ap_ufixed<3,3> dr2_msbs = dr2_full; // Convert via truncation
        ap_ufixed<2*POS_FRAC_BITS,0> dr2 = dr2_full;

        bool miss = (dr2_msbs!=0) || (dr2==0);

        //assert(dr2<0 && dr2<1);

        // General stuff
        direction_t dxn;
        distance_t dr;
        normalise_vector_16_v2(dxn, dr, dx, dr2);
        distance_t wr = 1 - dr;

        ////////////////////////////////////////////////////
        // Conservative
        auto conservative_strength = con_table[ (hb.type,ob.type) ];
        ap_ufixed<CONSERVATIVE_INT_BITS+COMPONENT_FRAC_BITS,CONSERVATIVE_INT_BITS> conservative_force = conservative_strength * wr;
        
        /////////////////////////////////////////////////////
        // Dissipative
        velocity_delta_t dv;
#pragma HLS UNROLL
        l1 : for(int d=0; d<3; d++){
            dv[d]=ob.v[d] - hb.v[d]; // Opposite sense, so that forces are all add (no sub)
        }
        auto sqrt_dissipative_strength = sqrt_diss_table[ (hb.type,ob.type) ];
        auto sqrt_gamma_p = round_to_frac<COMPONENT_FRAC_BITS>(sqrt_dissipative_strength * wr);

        auto gamma_p = round_to_frac<COMPONENT_FRAC_BITS>(sqrt_gamma_p * sqrt_gamma_p);

        auto dx_dot_dv = round_to_frac<COMPONENT_FRAC_BITS>(dxn[0] * dv[0] + dxn[1] * dv[1] + dxn[2] * dv[2]);

        auto dissipative_force = round_to_frac<COMPONENT_FRAC_BITS>(gamma_p * dx_dot_dv);


        ///////////////////////////////////////////////////
        // Stochastic
        // Assumed that scaled_hash_run is already pre-scaled for dt, e.g. sqrt(24/dt)
        auto u = hash_rng(round_seed, hb.id, ob.id);
        auto random_force = round_to_frac<COMPONENT_FRAC_BITS>(
            u* sqrt_gamma_p
        );

        //////////////////////////////////////////////////////
        // Final force
        auto total_force = conservative_force + dissipative_force + random_force;
        if(miss){
            total_force=0;
        }

#pragma HLS UNROLL
        for(int d=0; d<3; d++){
            hb_f[d] = round_to_frac<COMPONENT_FRAC_BITS>(
                dxn[d]*total_force
            );
        }
    }
};


void calc_force_fix_v2(
        const fixed_maths_v2::bead_info_t &hb,
        const fixed_maths_v2::bead_info_t &ob,
        ap_fixed_vec3<2,1> ob_cell_delta,
        const ap_uint<64> &round_seed,
        fixed_maths_v2::force_t &hb_f
) {
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS PIPELINE
    fixed_maths_v2::calc_force_fix(
        hb, ob, ob_cell_delta, round_seed,
        hb_f
    );

}


#endif
