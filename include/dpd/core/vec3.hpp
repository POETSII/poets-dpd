#ifndef vec3_hpp
#define vec3_hpp

// All good projects need yet another vector class...

#include <array>
#include <cmath>
#include <cassert>
#ifndef PDPD_TINSEL
#include <sstream>
#endif

#include "dpd_maths_primitives.hpp"

template<class T>
struct vec3g_t
{
    std::array<T,3> x{};

    vec3g_t()
        : x{0,0,0}
    {}

    vec3g_t(const vec3g_t &) = default;

    vec3g_t(T a, T b, T c)
        : x{a,b,c}
    {}

    template<class TT>
    explicit vec3g_t(const vec3g_t<TT> &_x)
        : x{(T)_x[0], (T)_x[1], (T)_x[2]}
    {}

    vec3g_t(const T _x[3])
        : x{_x[0], _x[1], _x[2]}
    {}

    template<class TT>
    void assign(const TT _x[3])
    {
        for(int i=0; i<3; i++){
            x[i]=_x[i];
        }
    }

    T operator[](size_t i) const
    { return x[i]; }

    T &operator[](size_t i)
    { return x[i]; }

    template<class TT>
    void extract(TT dst[3]) const
    {
        std::copy(x.begin(), x.end(), dst);
    }

    void clear()
    {
        x[0]=0;
        x[1]=0;
        x[2]=0;
    }

    template<class F>
    vec3g_t apply(T o, F f) const
    { return {f(x[0],o), f(x[1],o), f(x[2],o)}; }

    template<class F>
    vec3g_t apply(const vec3g_t &o, F f) const
    { return {f(x[0],o[0]), f(x[1],o[1]), f(x[2],o[2])}; }

    bool operator==(const vec3g_t &b) const
    { return x==b.x; }

    bool operator!=(const vec3g_t &b) const
    { return x!=b.x; }

    template<class O>
    vec3g_t operator+(const O &b) const
    { return apply(b, [](T aa, T bb){ return aa+bb; }); }

    template<class O>
    vec3g_t &operator+=(const O &x)
    { *this = *this + x; return *this; }

    template<class O>
    vec3g_t &operator-=(const O &x)
    { *this = *this - x; return *this; }

    template<class O>
    vec3g_t operator-(const O &b) const
    { return apply(b, [](T aa, T bb){ return aa-bb; }); }

    vec3g_t operator*(T b) const
    { return apply(b, [](T aa, T bb){ return aa*bb; }); }

    vec3g_t operator-() const
    { return { -x[0], -x[1], -x[2] }; }

    //! \note This is _not_ intended to be efficient
    vec3g_t reduce_to_box(const vec3g_t &box) const
    {
        vec3g_t res=*this;
        for(unsigned i=0; i<3; i++){
            while(res[i] < 0){
                res[i] += box[i];
            }
            while(res[i] >= box[i]){
                res[i] -= box[i];
            }
        }
        return res;
    }

    T dot(const vec3g_t &o) const
    {
        return x[0]*o.x[0] + x[1] *o.x[1] + x[2]*o.x[2];
    }

    T l2_norm() const
    {
        return pow_half(x[0]*x[0] + x[1]*x[1] + x[2] * x[2]);
    }
};

#ifndef PDPD_TINSEL
template<class T>
std::ostream &operator<<(std::ostream &dst, const vec3g_t<T> &x)
{
    return dst<<"("<<x[0]<<","<<x[1]<<","<<x[2]<<")";
}
#endif

using vec3r_t = vec3g_t<double>;
using vec3f_t = vec3g_t<float>;
using vec3i_t = vec3g_t<int>;

inline vec3r_t normalise(const vec3r_t &x)
{
    double r=x.l2_norm();
    if(r==0){
        return x;
    }else{
        return x * (1/r);
    }
}

inline bool isfinite(const vec3r_t &x)
{
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

inline bool isfinite(const vec3f_t &x)
{
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

inline double angle(const vec3r_t &a, const vec3r_t &b)
{
    double dot=a.dot(b);
    double l1=a.l2_norm();
    double l2=b.l2_norm();
    double cosPhi=dot / ( l1 * l2 );
    // It is fairly common for cosPhi to be slightly greater than 1 due to rounding.
    double res= (cosPhi >= 1) ? 0 : acos(cosPhi);
    assert(!std::isnan(res));
    return res;
}

inline vec3i_t floor(const vec3r_t &x)
{
    return {(int)std::floor(x[0]),(int)std::floor(x[1]),(int)std::floor(x[2])}; 
}

inline vec3i_t floor(const vec3f_t &x)
{
    return {(int)std::floor(x[0]),(int)std::floor(x[1]),(int)std::floor(x[2])}; 
}

/*
Wraps a vector into the range [0,bounds) in each dimension.
*/
template<class T>
inline vec3g_t<T> vec_wrap(const vec3g_t<T> &x, const vec3g_t<T> &bounds)
{
    return {
        x[0] + ((x[0]<0) ? bounds[0] : 0) + (x[0]>=bounds[0] ? -bounds[0] : 0),
        x[1] + ((x[1]<0) ? bounds[1] : 0) + (x[1]>=bounds[1] ? -bounds[1] : 0),
        x[2] + ((x[2]<0) ? bounds[2] : 0) + (x[2]>=bounds[2] ? -bounds[2] : 0)
    };
}

template<class T, class TC>
inline T vec3_wrapped_distance(const vec3g_t<T> &a, const vec3g_t<T> &b, const vec3g_t<TC> &box)
{
    T dist_sqr=0;
    for(unsigned i=0; i<3; i++){
        T dist=a[i]-b[i];
        if(dist<0){
            dist=-dist;
        }
        if(2*dist > box[i]){
            dist -= box[i];
        }
        dist_sqr += dist*dist;
    }
    return sqrt(dist_sqr);
}

namespace std
{
    template<>
    struct hash<vec3i_t>
    {
        size_t operator()(const vec3i_t &x) const
        {
            size_t res=5137385362063533605ull*(unsigned)x[0] + 17904276594807051537ull*(unsigned)x[1] + 17656613762115436855ull*(unsigned)x[2];
            return res;
        }
    };
}


template<class T>
inline void vec3_add(T x[3], const T y[3])
{ for(int i=0; i<3; i++){ x[i] += y[i]; } }

template<class T>
inline void vec3_add(T dst[3], const T x[3], const T y[3])
{ for(int i=0; i<3; i++){ dst[i] = x[i] + y[i]; } }

template<class T>
inline void vec3_sub(T dst[3], const T x[3])
{ for(int i=0; i<3; i++){ dst[i] -= x[i]; } }

template<class T>
inline void vec3_sub(T dst[3], const T x[3], const T y[3])
{ for(int i=0; i<3; i++){ dst[i] = x[i] - y[i]; } }

template<class T>
inline void vec3_mul(T dst[3], T s)
{ for(int i=0; i<3; i++){ dst[i] *= s; } }

template<class T>
inline void vec3_mul(T dst[3], const T x[3], T s)
{ for(int i=0; i<3; i++){ dst[i] = x[i] * s; } }

template<class T>
inline void vec3_neg(T dst[3])
{ for(int i=0; i<3; i++){ dst[i] = -dst[i]; } }

template<class T>
inline void vec3_add_mul(T dst[3], const T x[3], T s)
{ for(int i=0; i<3; i++){ dst[i] += x[i] * s; } }

template<class T>
inline T vec3_dot(const T x[3], const T y[3])
{ return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]; }

template<class T>
inline T vec3_dot_self(const T x[3])
{ return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; }

template<class T>
inline T vec3_l2_norm(const T x[3])
{
    return pow_half(vec3_dot_self(x));
}

template<class T>
inline void vec3_floor(int32_t dst[3], const T x[3])
{ for(int i=0; i<3; i++){ dst[i] = floor(x[i]); } }

// Does a floor with the pre-condition that x is non-negative
template<class T>
inline void vec3_floor_nn(int32_t dst[3], const T x[3])
{ for(int i=0; i<3; i++){ assert(x[i]>=0); dst[i] = floor_nn(x[i]); } }



template<class D, class T>
inline void vec3_copy(D &dst, const T &src)
{
    for(unsigned i=0; i<3; i++){
        dst[i]=src[i];
    }
}

template<class T>
inline bool vec3_equal(const T *x, const T *y)
{
    return x[0]==y[0] && x[1]==y[1] && x[2]==y[2];
}

template<class T>
inline void vec3_clear(T x[3])
{ for(int i=0; i<3; i++){ x[i]=0; } }

template<class T>
inline bool vec3_is_zero(const T *x)
{ return x[0]==0 && x[1]==0 && x[2]==0; }


#endif
