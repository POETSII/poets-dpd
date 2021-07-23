#ifndef test_multiple_moving_beads_hpp
#define test_multiple_moving_beads_hpp

#include "test_base.hpp"

#include "../dpd_state_builder.hpp"

#include <iostream>

class TestMultipleMovingBeads
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    std::string m_name;
    std::vector<vec3r_t> m_x0;
    vec3r_t m_v;
    vec3r_t m_sep;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx, const vec3r_t &sep)
    {
        std::stringstream tmp;
        tmp<<"TestMultipleMovingBeads[x0="<<x0<<";dx="<<dx<<";sep="<<sep<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [](vec3r_t dx)
        {
            vec3r_t sep{1.2, 1.2, 1.2};

            std::string name=calc_name({0,0,0}, dx, sep);
            TestBase::add_test_factory(name, 0, [=](){
                return std::make_shared<TestMultipleMovingBeads>(name,  dx, sep);
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

    TestMultipleMovingBeads(std::string name, const vec3r_t &v=vec3r_t{0.5,0,0}, const vec3r_t &sep=vec3r_t{1.2, 1.2, 1.2})
        : m_name(name)
        , m_v(v)
        , m_sep(sep)
    {}

    virtual std::string name() const 
    {
        return m_name;
    }

    WorldState create_world() override
    {
        vec3r_t len={7,8,9};

        WorldStateBuilder b(len);
        WorldState &s=b.data();

        s.t=0;
        s.dt=1.0/64;

        b.add_bead_type("W");
        b.add_polymer_type("W", {"W"}, {}, {});

        vec3r_t pos={0,0,0};
        while(pos[0] + m_sep[0] < len[0]){
            while(pos[1] + m_sep[1] < len[1]){
                while(pos[2] + m_sep[2] < len[2]){
                    unsigned index=b.add_monomer("W", pos, m_v);
                    assert(index==m_x0.size());
                    m_x0.push_back(pos);
                    pos[2] += m_sep[2];
                }
                pos[2]=0;
                pos[1] += m_sep[1];
            }
            pos[1]=0;
            pos[0] += m_sep[0];
        }

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        
        if( abs(s.t - s.dt * m_steps_done ) > 1e-10 ){
            throw TestFailedException("t is not dt*nSteps.");
        }
        for(unsigned i=0; i<m_x0.size(); i++){
            vec3r_t pos=s.beads.at(i).x;
            vec3r_t full_x=m_v * s.t + m_x0[i];
            vec3r_t red_x=full_x.reduce_to_box(s.box);
            //std::cerr<<"Pos = "<<s.beads.at(0).x<<"\n";
            double err= (red_x - pos ).l2_norm();
            if(err > sqrt(m_steps_done) * 1e-4){
                std::stringstream acc;
                acc<<"Bead is in the wrong position at t="<<s.t<<". Expected "<<red_x<<", got "<<pos<<", err="<<err<<", tol="<<sqrt(m_steps_done) * 1e-4;
                throw TestFailedException(acc.str());
            }
        }

        if(s.t < 10){
            m_steps_done += m_step_dist;
            return m_step_dist++;
        }else{
            return 0;
        }
    }
};

#endif
