#include "hash.hpp"

#include <string>
#include <memory>
#include <cmath>
#include <cassert>

extern "C"{
    #include <testu01/bbattery.h>
    #include <testu01/swrite.h>
}

struct enum_fixed_pair_step_t
{
    uint32_t a;
    uint32_t b;
    uint64_t t;

    std::string name() const
    {
        return "FixedPairStepT[a="+std::to_string(a)+";b="+std::to_string(b)+"]";
    }

    uint32_t operator()()
    {
        auto t_hash=next_t_hash(t);
        return hash_rng_sym(t_hash, a, b);
    }
};

struct enum_fixed_a_t_inc_b
{
    uint32_t a;
    uint32_t b;
    uint64_t t_hash;

    enum_fixed_a_t_inc_b(uint32_t _a, uint32_t _b, uint64_t t)
        : a(_a)
        , b(_b)
    {
        t_hash = next_t_hash(t);
    }

    std::string name() const
    { return "FixedATIncB[a="+std::to_string(a)+";t="+std::to_string(t_hash)+"]";
    }

    uint32_t operator()()
    {
        b++;
        return hash_rng_sym(t_hash, a, b);
    }
};

struct enum_fixed_t_inc_a_b
{
    uint32_t a;
    uint32_t b;
    uint64_t t_hash;

    enum_fixed_t_inc_a_b(uint32_t _a, uint32_t _b, uint64_t t)
        : a(_a)
        , b(_b)
    {
        t_hash = next_t_hash(t);
    }


    std::string name() const
    { return "FixedTIncAB[t="+std::to_string(t_hash)+"]";
    }

    uint32_t operator()()
    {
        a++;
        b++;
        return hash_rng_sym(t_hash, a, b);
    }
};

struct enum_wrap
{
    uint32_t n;
    uint32_t a;
    uint32_t b;
    uint64_t t_hash;
    uint64_t t;

    enum_wrap(uint32_t _n)
        : n(_n)
        , a(0)
        , b(0)
        , t(123456789)
    {
        t_hash=next_t_hash(t);
    }

    std::string name() const
    { return "1K";
    }

    uint32_t operator()()
    {
        ++b;
        if(a<=b){
            b=0;
            ++a;
            if(a==n){
                a=1;
                t_hash=next_t_hash(t);
            }
        }
        assert(b<a);
        //fprintf(stderr, "%08x, %08x, %08x\n", t_hash, a, b);
        return hash_rng_sym(t_hash, a, b);
    }
};

struct enum_wrap_rev
{
    uint32_t n;
    uint32_t a;
    uint32_t b;
    uint64_t t_hash;
    uint64_t t;

    enum_wrap_rev(uint32_t _n)
        : n(_n)
        , a(0)
        , b(0)
        , t(123456789)
    {
        t_hash=next_t_hash(t);
    }

    std::string name() const
    { return "1Krev";
    }

    uint32_t rev(uint32_t v)
    {
        // http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv
        // swap odd and even bits
        v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
        // swap consecutive pairs
        v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
        // swap nibbles ... 
        v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
        // swap bytes
        v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
        // swap 2-byte long pairs
        v = ( v >> 16             ) | ( v               << 16);
        return v;
    }

    uint32_t operator()()
    {
        ++b;
        if(a<=b){
            b=0;
            ++a;
            if(a==n){
                a=1;
                t_hash=next_t_hash(t);
            }
        }
        assert(b<a);
        //fprintf(stderr, "%08x, %08x, %08x\n", t_hash, a, b);
        return hash_rng_sym(t_hash, rev(a), rev(b));
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

    auto gen=make_gen( enum_fixed_pair_step_t{1,1,1} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_fixed_a_t_inc_b{1,1,1} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_fixed_t_inc_a_b{1,1,2} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap_rev{10} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap{10} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap{100} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap{1000} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap{10000} );
    bbattery_SmallCrush(gen.get());

    gen=make_gen( enum_wrap{100000} );
    bbattery_SmallCrush(gen.get());
}