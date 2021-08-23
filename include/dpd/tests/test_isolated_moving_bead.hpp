#ifndef test_isolated_moving_bead_hpp
#define test_isolated_moving_bead_hpp

#include "test_base.hpp"

#include <iostream>
#include <random>

class TestIsolatedMovingBead
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    std::string m_name;
    vec3r_t m_v;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx)
    {
        std::stringstream tmp;
        tmp<<"TestIsolatedMovingBead[x0="<<x0<<";dx="<<dx<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [](vec3r_t dx)
        {
            std::string name=calc_name({0,0,0}, dx);
            TestBase::add_test_factory(name, 0, [=](){
                return std::make_shared<TestIsolatedMovingBead>(name,  dx);
            });
        };

        vec3r_t v;
        for(v[0]=-1; v[0]<=1; v[0]++){
            for(v[1]=-1; v[1]<=1; v[1]++){
                for(v[2]=-1; v[2]<=1; v[2]++){
                    add(v);
                }
            }
        }

        std::mt19937_64 rng(1);
        std::uniform_real_distribution<> urng(-1,1);
        for(int i=0; i<10; i++){
            add({
                urng(rng), urng(rng), urng(rng)
            });
        }

    }

    TestIsolatedMovingBead(std::string name, const vec3r_t &v={0.5,0,0})
        : m_name(name)
        , m_v(v)
    {}

    virtual std::string name() const 
    {
        return m_name;
    }

    WorldState create_world() override
    {
        WorldState res;
        for(unsigned i=0; i<3; i++){
            res.origin[i]=0;
            res.box[i]=4;
        }
        res.t=0;
        res.dt=1.0/64;

        res.interactions={
            {1, 1}
        };
        res.bead_types.push_back({
            "W", 0.5, 0
        });
        res.polymer_types.push_back({
            "WW",
            { 0 },
            {},
            {},
            0
        });
        res.beads.push_back({
            0, 0, 0,
            0, 0, true,
            { 0, 0, 0 },
            {m_v.x[0], m_v.x[1], m_v.x[2]},
            { 0, 0, 0 }
        });
        res.polymers.push_back({
            {0}, 0, 0
        });

        return res;
    }

    unsigned get_advance_count(WorldState &s) override
    {
        vec3r_t pos=s.beads.at(0).x;
        vec3r_t full_x=m_v * s.t * s.dt;
        vec3r_t red_x=full_x.reduce_to_box(s.box);
        //std::cerr<<"Pos = "<<s.beads.at(0).x<<"\n";
        double err= (red_x - pos ).l2_norm();
        if(err > sqrt(m_steps_done) * 1e-5){
            std::stringstream acc;
            acc<<"Bead is in the wrong position at t="<<s.t<<". Expected "<<red_x<<", got "<<pos;
            throw TestFailedException(acc.str());
        }

        if(s.t*s.dt < 10){
            m_steps_done += m_step_dist;
            return m_step_dist++;
        }else{
            return 0;
        }
    }
};

#endif
