#ifndef fixed_math_hpp
#define fixed_math_hpp

#include "ap_fixed_vec3.hpp"

ap_uint<32> xorshift(ap_uint<32> x)
    {
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return x;
    };

ap_fixed<17,1> hash_rng(ap_uint<64> round_seed, ap_uint<32> hb_id, ap_uint<32> ob_id)
{
    auto lt = hb_id < ob_id;
    auto a = lt ? hb_id : ob_id;
    auto b = lt ? ob_id : hb_id;

    auto left=xorshift( round_seed.range(63,32) ^ a);
    auto right=xorshift( round_seed.range(31,0) ^ b); // Should be a different xorshift
    auto mid=xorshift( left ^ right );
    mid=xorshift(mid);

    // TODO : this should be scaled by sqrt(24/dt)
    return reinterpret_as_ap_fixed<17,1>(mid);
}

struct fixed_maths
{
    static const int BEAD_ID_BITS = 32;
    static const int ROUND_SEED_BITS = 32;

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

    __attribute__((always_inline)) static void calc_force_fix(
        const ap_uint<BEAD_ID_BITS> &hb_id,
        const ap_uint<BEAD_ID_BITS> &ob_id,
        const ap_uint<ROUND_SEED_BITS> &round_seed,
        const conservative_strength_t &conservative_strength,
        const sqrt_dissipative_strength_t &sqrt_dissipative_strength,
        const pos_delta_t &dx,  // Each component strictly in (-1,+1)
        const distance_sqr_t &dr2,   // Strictly in (0,1)  (Cannot be 0 as would be self-interaction)
        const velocity_t &hb_v,
        const velocity_t &ob_v,
        force_t &hb_f
    ){
        //assert(dr2<0 && dr2<1);

        // General stuff
        direction_t dxn;
        distance_t dr;
        normalise_vector_16_v2(dxn, dr, dx, dr2);
        distance_t wr = 1 - dr;

        ////////////////////////////////////////////////////
        // Conservative
        ap_ufixed<CONSERVATIVE_INT_BITS+COMPONENT_FRAC_BITS,CONSERVATIVE_INT_BITS> conservative_force = conservative_strength * wr;
        
        /////////////////////////////////////////////////////
        // Dissipative
        velocity_delta_t dv;
#pragma HLS UNROLL
        l1 : for(int d=0; d<3; d++){
            dv[d]=ob_v[d] - hb_v[d]; // Opposite sense, so that forces are all add (no sub)
        }
        auto sqrt_gamma_p = round_to_frac<COMPONENT_FRAC_BITS>(sqrt_dissipative_strength * wr);

        auto gamma_p = round_to_frac<COMPONENT_FRAC_BITS>(sqrt_gamma_p * sqrt_gamma_p);

        auto dx_dot_dv = round_to_frac<COMPONENT_FRAC_BITS>(dxn[0] * dv[0] + dxn[1] * dv[1] + dxn[2] * dv[2]);

        auto dissipative_force = round_to_frac<COMPONENT_FRAC_BITS>(gamma_p * dx_dot_dv);


        ///////////////////////////////////////////////////
        // Stochastic
        // Assumed that scaled_hash_run is already pre-scaled for dt, e.g. sqrt(24/dt)
        auto u = hash_rng(round_seed, hb_id, ob_id);
        auto random_force = round_to_frac<COMPONENT_FRAC_BITS>(
            u* sqrt_gamma_p
        );

        //////////////////////////////////////////////////////
        // Final force
        auto total_force = conservative_force + dissipative_force + random_force;

#pragma HLS UNROLL
        for(int d=0; d<3; d++){
            hb_f[d] = round_to_frac<COMPONENT_FRAC_BITS>(
                dxn[d]*total_force
            );
        }
    }
};


void calc_force_fix(
        const ap_uint<32> &hb_id,
        const ap_uint<32> &ob_id,
        const ap_uint<64> &round_seed,
        const fixed_maths::conservative_strength_t &conservative_strength,
        const fixed_maths::sqrt_dissipative_strength_t &sqrt_dissipative_strength,
        const ap_fixed<17,1> &dx_0,  // Each component strictly in (-1,+1)
		const ap_fixed<17,1> &dx_1,  // Each component strictly in (-1,+1)
		const ap_fixed<17,1> &dx_2,  // Each component strictly in (-1,+1)
		const fixed_maths::distance_sqr_t &dr2,   // Strictly in (0,1)  (Cannot be 0 as would be self-interaction)
        const fixed_maths::velocity_part_t &hb_v_0,
		const fixed_maths::velocity_part_t &hb_v_1,
		const fixed_maths::velocity_part_t &hb_v_2,
		const fixed_maths::velocity_part_t &ob_v_0,
		const fixed_maths::velocity_part_t &ob_v_1,
		const fixed_maths::velocity_part_t &ob_v_2,
        fixed_maths::force_part_t &hb_f_0,
		fixed_maths::force_part_t &hb_f_1,
		fixed_maths::force_part_t &hb_f_2
)
{
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS PIPELINE
	fixed_maths::force_t hb_f;
    fixed_maths::calc_force_fix(
        hb_id, ob_id, round_seed,
        conservative_strength, sqrt_dissipative_strength,
        {dx_0,dx_1,dx_2}, dr2,
        {hb_v_0,hb_v_1,hb_v_2}, {ob_v_0,ob_v_1,ob_v_2},
        hb_f
    );
    hb_f_0=hb_f[0];
    hb_f_1=hb_f[1];
    hb_f_2=hb_f[2];

}


#endif
