#ifndef ap_fixed_vec3_hpp
#define ap_fixed_vec3_hpp

#include <ap_fixed.h>
#include <hls_math.h>
#include <array>

template<int MSB,int LSB,int Width,int IntBits,ap_q_mode QMode>
ap_uint<MSB-LSB+1> get_bits(ap_ufixed<Width,IntBits,QMode> x)
{
    return x.range(MSB,LSB);
}

template<int Width, int IntBits, ap_q_mode QMode=AP_TRN>
using ap_fixed_vec3 = std::array<ap_fixed<Width,IntBits,QMode>,3>;

template<int Width, int IntBits, ap_q_mode QMode=AP_TRN>
ap_fixed<Width,IntBits,QMode> reinterpret_as_ap_fixed(ap_uint<Width> x)
{
    ap_fixed<Width,IntBits,QMode> res;
    res.V=x;
    return res;
}

template<int FracBits,int Width,int IntBits,ap_q_mode QMode>
ap_fixed<IntBits+FracBits,IntBits,QMode> round_to_frac(ap_fixed<Width,IntBits,QMode> x)
{
    //static_assert( (Width-IntBits) > FracBits);
    auto t1=ap_fixed<IntBits+FracBits,IntBits,AP_RND_CONV>(x);
    return ap_fixed<IntBits+FracBits,IntBits,QMode>(t1);
}

template<int FracBits,int Width,int IntBits,ap_q_mode QMode>
ap_ufixed<IntBits+FracBits,IntBits,QMode> round_to_frac(ap_ufixed<Width,IntBits,QMode> x)
{
    //static_assert( (Width-IntBits) > FracBits);
    auto t1=ap_ufixed<IntBits+FracBits,IntBits,AP_RND_CONV>(x);
    return ap_ufixed<IntBits+FracBits,IntBits,QMode>(t1);
}

template<int FracBits,int Width,int IntBits,ap_q_mode QMode>
ap_ufixed<IntBits+FracBits,IntBits,QMode> trunc_to_frac(ap_ufixed<Width,IntBits,QMode> x)
{
    //static_assert( (Width-IntBits) > FracBits);
    auto t1=ap_ufixed<IntBits+FracBits,IntBits,AP_TRN>(x);
    return ap_ufixed<IntBits+FracBits,IntBits,QMode>(t1);
}


// Pre: dr2 in (0,1)
void normalise_vector_16(
    ap_fixed_vec3<1+16,1> &dxn, ap_ufixed<16,0> &r,
    ap_fixed_vec3<1+16,1> dx, ap_ufixed<16*2,0> dr2
){
    // dr2 is in [2^-32,1), 1/dr2 is in [2^16,1)
    //ap_ufixed<24,17,AP_RND_CONV> inv_sqrt_dr2=hls::rsqrt(ap_ufixed<24,17>(dr2));
    // ERROR : This is now complete nonsense
	ap_ufixed<24,17,AP_RND_CONV> inv_sqrt_dr2 = dr2;
	inv_sqrt_dr2 = inv_sqrt_dr2 + 7 * inv_sqrt_dr2 * inv_sqrt_dr2;
	inv_sqrt_dr2 = inv_sqrt_dr2 + 7 * inv_sqrt_dr2 * inv_sqrt_dr2;

    for(int i=0;i<3;i++){
        dxn[i] = inv_sqrt_dr2 * ap_fixed<1+16,1,AP_RND_CONV>(dx[i]);
    }
    r = inv_sqrt_dr2 * ap_ufixed<32,0,AP_RND_CONV>(dr2);
}


// Pre: dr2 in (0,1)
void normalise_vector_16_v2(
    ap_fixed_vec3<1+16,1> &dxn, ap_ufixed<16,0> &r,
    ap_fixed_vec3<1+16,1> dx, ap_ufixed<16*2,0> dr2
){
    // dr2 is in [2^-32,1), 1/dr2 is in [2^16,1)
    ap_ufixed<16*2,0> dr2n=dr2;
    if(get_bits<31,16>(dr2)==0){
        dr2 <<= 16;
        for(int d=0; d<3; d++){
            dx[d] <<= 8;
        }
        dr2n <<= 8;
    }
    if(get_bits<31,24>(dr2)==0){
        dr2 <<= 8;
        for(int d=0; d<3; d++){
            dx[d] <<= 4;
        }
        dr2n <<= 4;
    }
    if(get_bits<31,28>(dr2)==0){
        dr2 <<= 4;
        for(int d=0; d<3; d++){
            dx[d] <<= 2;
        }
        dr2n <<= 2;
    }
    if(get_bits<31,30>(dr2)==0){
        dr2 <<= 2;
        for(int d=0; d<3; d++){
            dx[d] <<= 1;
        }
        dr2n <<= 1;
    }

    // x in the range [0.25,1)
    // y in the range [4,1), so 
    auto half_x = ap_ufixed<18,-1>(dr2>>1); 
    double exp=1.0/sqrt(dr2.to_double());

    // Fake table with 64 entries
    //auto index=ap_ufixed<6,0,AP_TRN>(dr2);
    //double approx = 1.0/sqrt(index.to_double());
    auto index=get_bits<31,26>(dr2);
    static const ap_ufixed<16,4> table[64]={
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.969463856,1.912365775,1.85996222,1.811643255,1.766904417,1.725324371,1.686548085,1.650273994,1.616244071,1.584236069,1.55405738,1.525540143,1.498537299,1.472919389,1.448571937,1.42539329,1.403292831,1.382189481,1.362010449,1.342690173,1.324169422,1.306394529,1.289316742,1.272891655,1.257078722,1.241840841,1.227143982,1.21295687,1.199250702,1.185998907,1.17317692,1.160762,1.148733054,1.137070487,1.125756072,1.114772823,1.104104895,1.093737483,1.083656738,1.073849688,1.064304168,1.055008757,1.045952721,1.037125958,1.028518954,1.020122741,1.011928851,1.003929288
    };
    ap_ufixed<12,2> y0 = table[index];
    //std::cerr<<"  dr2 = "<<dr2<<"\n";
    //std::cerr<<"  exp="<<exp<<", y0="<<y0<<"\n";
    ap_ufixed<17,2> y1 = y0 * (ap_ufixed<2,1>(1.5) - trunc_to_frac<15>(trunc_to_frac<15>(half_x * y0) * y0)  );
    //std::cerr<<"  exp="<<exp<<", y1="<<y1<<"\n";
    ap_ufixed<20,2> y2 = y1 * (ap_ufixed<2,1>(1.5) - round_to_frac<18>(round_to_frac<18>(half_x * y1) * y1)  );
    //std::cerr<<"  exp="<<exp<<", y2="<<y2<<"\n";
    
    for(int i=0;i<3;i++){
        dxn[i] = ap_fixed<1+16,1,AP_RND_CONV,AP_SAT>(y2 * ap_fixed<1+16,1,AP_RND_CONV>(dx[i]));
    }
    r = ap_ufixed<16,0,AP_RND_CONV,AP_SAT>( y2 * ap_ufixed<18,0,AP_RND_CONV,AP_SAT>(dr2n) );
}

#endif
