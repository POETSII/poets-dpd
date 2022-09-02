#include "dpd/core/hash.hpp"

#include <string>
#include <memory>
#include <cmath>
#include <cassert>
#include <random>

extern "C"{
    #include <testu01/bbattery.h>
    #include <testu01/swrite.h>
}

__m256i f(__m256i &x)
{
    __m256i a=x;

    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
    x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));

    __m256i b=_mm256_shuffle_epi32(x, (2<<6)|(3<<4)|(0<<2)|(1<<0));

    /*x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
    x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));    

    __m256i c=x;*/

    return _mm256_add_epi32(a,b);
}

__m256i f2(__m256i &x, __m256i &r)
{
    __m256i a=x;

    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
    x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));

    __m256i b=x;

    r = _mm256_add_epi32(r, _mm256_set1_epi32(19937));  

    return _mm256_add_epi32(x, r);
}

struct enum_horizontal_t
{
    __m256i state;
    uint32_t values[8];
    unsigned n;

    enum_horizontal_t(std::mt19937_64 &u)
    {
        for(unsigned i=0; i<4; i++){
            state[i] = u();
        }
        n=0;
    }

    std::string name() const
    {
        return "Horizontal";
    }

    uint32_t operator()()
    {
        if(n==0){
            __m256i res=f(state);
            _mm256_storeu_si256((__m256i*)values, res);
            n=8;
        }
        return values[--n];
    }
};

struct enum_vertical_t
{
    __m256i state;
    uint32_t values[8];
    unsigned index;

    enum_vertical_t(std::mt19937_64 &u, unsigned sel=0)
    {
        for(unsigned i=0; i<4; i++){
            state[i] = u();
        }
        index=sel;
    }

    std::string name() const
    {
        return "Vertical";
    }

    uint32_t operator()()
    {
        __m256i res=f(state);
        _mm256_storeu_si256((__m256i*)values, res);
        return values[index];
    }
};

struct enum_horizontal_f2_t
{
    __m256i state[2];
    uint32_t values[8];
    unsigned n;

    enum_horizontal_f2_t(std::mt19937_64 &u)
    {
        for(unsigned i=0; i<4; i++){
            state[0][i] = u();
            state[1][i] = u();
        }
        n=0;
    }

    std::string name() const
    {
        return "HorizontalF2";
    }

    uint32_t operator()()
    {
        if(n==0){
            __m256i res=f2(state[0], state[1]);
            _mm256_storeu_si256((__m256i*)values, res);
            n=8;
        }
        return values[--n];
    }
};

template<class T>
struct adapter
    : unif01_Gen
{
    T enumerator;
    std::string sname;

    static unsigned long next_bits(void *, void *s)
    {
        return ((adapter*)s)->enumerator();
    }

    static double next_u01(void *, void *s)
    {
        const double scale=ldexp(1, -32);
        return ((adapter*)s)->enumerator() * scale;
    }

    static void write(void *s)
    {
        fprintf(stdout, "%s", ((adapter*)s)->name);
    }

    adapter(const T &e)
        : enumerator(e)
        , sname(e.name())
    {
        this->GetBits=next_bits;
        this->GetU01=next_u01;
        this->name=(char*)sname.c_str();
        this->param=0;
        this->state=this;
        this->Write=write;
    }
};


template<class E>
std::unique_ptr<unif01_Gen> make_gen(const E &x)
{
    return std::make_unique<adapter<E>>(x);
}

int main()
{
    swrite_Basic=true;

    std::mt19937_64 init;

    auto gen=make_gen( enum_horizontal_t(init) );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_vertical_t(init, 0) );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_vertical_t(init, 1) );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_horizontal_t(init) );
    bbattery_Crush(gen.get());

    gen=make_gen( enum_vertical_t(init, 0) );
    bbattery_Crush(gen.get());

    gen=make_gen( enum_vertical_t(init, 1) );
    bbattery_Crush(gen.get());
}