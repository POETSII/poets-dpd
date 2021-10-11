#ifndef CVec3_hpp
#define CVec3_hpp

#include "vec3.hpp"

#include <cmath>

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
            m_exp=0;
            for(int i=0; i<3; i++){
                m_p[i]=0;
            }
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
            frs[i]=round( ldexp(fs[i], -e+15 ) );
            overflow=overflow||(std::abs(frs[i])==(1<<15));
        }

        if( overflow ){
            // One component got rounded up to 1.0
            e=e+1;
            for(int i=0; i<3; i++){
                frs[i]=round( ldexp(fs[i], -e+15 ) );
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
        double scale=ldexp(1.0/32768, m_exp);
        return {m_p[0]*scale, m_p[1]*scale, m_p[2]*scale};
    }
};

static_assert(sizeof(CVec3)==8);

#endif
