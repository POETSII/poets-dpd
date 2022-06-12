#include "dpd/core/hash.hpp"

#include <string>
#include <memory>
#include <cmath>
#include <cassert>

extern "C"{
    #include <testu01/bbattery.h>
    #include <testu01/swrite.h>
}

__m256i f(__m256i &x)
{
    __m256i ox=x;

    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 13));
    x=_mm256_xor_si256(x, _mm256_srli_epi64(x, 7));
    x=_mm256_xor_si256(x, _mm256_slli_epi64(x, 17));    

    auto y=_mm256_shuffle_epi32(x, (2<<6)|(3<<4)|(0<<2)|(1<<1));

    return _mm256_add_epi32(x,y);
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
        n=;
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
        return state[--n];
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
    swrite_Basic=false;

    std::mt19937_64 init;

    auto gen=make_gen( init );
    bbattery_SmallCrush(gen.get(enum_horizontal_t));
}