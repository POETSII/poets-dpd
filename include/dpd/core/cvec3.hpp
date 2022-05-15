#ifndef CVec3_hpp
#define CVec3_hpp

#include "vec3.hpp"

#include <cmath>
#include <iostream>

struct CVec3
{
private:
    // Each component is a fraction in the range (-1,1)
    // scaled by 2^15. So the raw values are (-2^15..2^15),
    // or equivalently [-2^15+1..2^15-1]
    //
    // If the number is non-zero then it is guaranteed that
    // at least one component is in the range (-1,-0.5] or [0.5,1)
    //
    // The exponent applies to each component equally. So given
    // a component c in (-1,1), the value is c*2^exp.
    // Or given a raw component r in (-2^15..2^15) the value is c*2^-15*2^exp==c*2^(exp-15)
    //
    // Zero is recorded as all components zero, and exp also at zero.
    // If any component is non-zero, the vector is non-zero.
    //
    // Inf and nan are not supported.

    int16_t m_exp;
    int16_t m_p[3];

    void set(const double po[3])
    {
        bool is_zero=true;
        double p[3];
        for(int i=0; i<3; i++){
            auto c=std::fpclassify(po[i]);
            switch(c){
            case FP_INFINITE: throw std::runtime_error("Attempt to convert nan or inf to CVec3.");
            case FP_NAN: throw std::runtime_error("Attempt to convert nan or inf to CVec3.");
            case FP_SUBNORMAL:
            case FP_ZERO: p[i]=0; break;
            case FP_NORMAL: p[i]=po[i]; is_zero=false; break;
            default: assert(0); break;
            }
        }

        if(is_zero){
            memset(this, 0, sizeof(CVec3));
            return;
        }


        int es[3];
        double fs[3];
        int e=INT16_MIN;
        for(int i=0; i<3; i++){
            fs[i]=frexp(p[i], &es[i]);
            e=std::max(e, es[i]);
        }

        double frs[3];
        bool overflow=false;
        for(int i=0; i<3; i++){
            frs[i]=round( ldexp(p[i], 15-e ) );
            overflow=overflow||(std::abs(frs[i])==(1<<15));
        }

        if( overflow ){
            // One component got rounded up to 1.0
            e=e+1;
            for(int i=0; i<3; i++){
                frs[i]=round( ldexp(p[i], 15-e ) );
            }
        }
        
        m_exp=e;
        for(int i=0; i<3;i++){
            m_p[i]=(int16_t)frs[i];
        }
    }

    double to_double(int16_t v) const
    { return ldexp( v, m_exp-15 ); }
public:
    CVec3()
    {
        memset(this, 0, sizeof(CVec3));
    }

    CVec3(double x, double y, double z)
    {
        double p[3]={x,y,z};
        set(p);
    }

    CVec3(const vec3r_t &v)
    {
        set(&v.x[0]);
    }


    double get_x() const
    { return to_double(m_p[0]); }

    double get_y() const
    { return to_double(m_p[1]); }

    double get_z() const
    { return to_double(m_p[2]); }

    vec3r_t get_vec3r() const
    {
        double scale=ldexp(1.0, m_exp-15);
        return {m_p[0]*scale, m_p[1]*scale, m_p[2]*scale};
    }
};

static_assert(sizeof(CVec3)==8);

struct CVec3Half
{
private:
    // Packs a direction vector into 4 bytes
    uint32_t exponent : 5; // Exponent in the range 0=0, 1=2^-14, 2=2^-12, ..., 30=2^15, 31=invalid
    int32_t fx : 9; // Value in (-1,+1) encoded as twos-complement integer with 8 fractional bits and a sign bit
    int32_t fy : 9;
    int32_t fz : 9;
    // 0 is represent by raw exponent field of 0. All other fields must be 0.

    void set(const double po[3])
    {
        bool is_zero=true;
        double p[3];
        for(int i=0; i<3; i++){
            auto c=std::fpclassify(po[i]);
            switch(c){
            case FP_INFINITE: throw std::runtime_error("Attempt to convert nan or inf to CVec3Half.");
            case FP_NAN: throw std::runtime_error("Attempt to convert nan or inf to CVec3Half.");
            case FP_SUBNORMAL:
            case FP_ZERO: p[i]=0; break;
            case FP_NORMAL: p[i]=po[i]; is_zero=false; break;
            default: assert(0); break;
            }
        }

        if(is_zero){
            memset(this, 0, sizeof(CVec3Half));
            return;
        }


        int es[3];
        double fs[3];
        int e=INT16_MIN;
        for(int i=0; i<3; i++){
            fs[i]=frexp(p[i], &es[i]);
            e=std::max(e, es[i]);
        }
        assert(e!=INT16_MIN); // Onecomponent is non-zero.

        if(e > 15){
            std::cerr<<"x="<<p[0]<<", y="<<p[1]<<", z="<<p[2]<<"\n";
            throw std::runtime_error("Vector overflowed while converting to CVec3Half, so magnitude is 65536 or more.");
        }
        if(e < -14){
            memset(this, 0, sizeof(CVec3Half));
            return;
        }

        double frs[3];
        bool overflow=false;
        for(int i=0; i<3; i++){
            frs[i]=round( ldexp(p[i], 8-e ) );
            overflow=overflow||(std::abs(frs[i])==(1<<8));
        }

        if( overflow ){
            // One component got rounded up to 1.0
            e=e+1;
            if(e > 15){
                std::cerr<<"x="<<p[0]<<", y="<<p[1]<<", z="<<p[2]<<"\n";
                throw std::runtime_error("Vector overflowed while converting to CVec3Half, so magnitude is 65536 or more.");
            }
            for(int i=0; i<3; i++){
                frs[i]=round( ldexp(p[i], 8-e ) );
            }
        }        
        exponent=e+15;
        fx=frs[0];
        fy=frs[1];
        fz=frs[2];

        //fprintf(stderr, "e=%d (=2^%d), fs[0]=%g, frs[0]=%g, fx=%d=(%f), fy=%d, fz=%d\n", exponent, exponent-15, fs[0], frs[0], fx, fx/256.0, fy, fz);
    }

    double to_double(int16_t v) const
    { return ldexp( v, exponent-15-8 ); }
public:
    CVec3Half()
    {
        memset(this, 0, sizeof(CVec3Half));
    }

    CVec3Half(double x, double y, double z)
    {
        double p[3]={x,y,z};
        set(p);
    }

    CVec3Half(const vec3r_t &v)
    {
        set(&v.x[0]);
    }


    double get_x() const
    { return to_double(fx); }

    double get_y() const
    { return to_double(fy); }

    double get_z() const
    { return to_double(fz); }

    vec3r_t get_vec3r() const
    {
        double scale=ldexp(1.0, exponent-15-8);
        return {fx*scale, fy*scale, fz*scale};
    }

    bool is_zero() const
    {
        assert( (exponent!=0) == ( (fx||fy||fz) ) );
        return exponent==0;
    }
};

static_assert(sizeof(CVec3Half)==4);

/*
struct CVec3Polar
{
private:
    // This is a compressed version of CVec3 that packs
    // into 32 bits using polar co-ordinates
    union{
        struct {
            uint32_t r_exp  :  6; // exponent, biased by 32
            uint32_t r_frac :  9; // fraction, in [0.5,1) with 1 implicit bit
            uint32_t phi    :  9; // angle
            uint32_t theta  :  8; // azimuth
        };
        uint32_t bits;
    };

    // Zero is recorded as all components zero.
    //
    // Inf and nan are not supported.

    static std::pair<double,double> sincos(uint angle)
    {
        assert(angle < 512);
        bool table_inited=false;
        static std::array<std::pair<double,double>,512> table;

        return table.at(angle);
    }
public:
    CVec3Polar()
    {
        bits=0;
    }

    CVec3Polar(double x, double y, double z)
    {
        double p[3]={x,y,z};
        set(p);
    }

    CVec3Polar(const vec3r_t &v)
    {
        set(&v.x[0]);
    }

    void set(const double p[3])
    {
        if(p[0]==0 && p[1]==0 && p[2]==0){
            bits = 0;
            return;
        }

        float full_r = sqrtf(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        float full_theta = atan2f(sqrt(p[0]*p[0] + p[1]*p[1]), p[2]);
        float full_phi = atan2f(y,x);

        int full_r_exp;
        full_r_frac = frexp( full_r, full_r_exp );

        TODO
    }


    vec3r_t get_vec3r() const
    {
        if(r_frac==0){
            return {0,0,0};
        }

        auto sc_theta=sincos(theta);
        auto sc_phi=sincos(phi);
        double full_frac = r_frac | (1u<<9);
        double r = ldexp(full_frac, r_exp-32-9 );

        return { r * sc_phi.second * sc_theta.first, r * sc_phi.first * sc_theta.first, r * sc_theta.second };
    }
};

static_assert(sizeof(CVec3Polar)==4);
*/

#endif
