#ifndef test_multiple_moving_beads_hpp
#define test_multiple_moving_beads_hpp

#include "test_base.hpp"

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
    vec3r_t m_dims;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx, const vec3r_t &sep, const vec3r_t &dims)
    {
        std::stringstream tmp;
        tmp<<"TestMultipleMovingBeads[x0="<<x0<<";dx="<<dx<<";sep="<<sep<<";dims="<<dims<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [](vec3r_t dx, vec3r_t dims)
        {
            vec3r_t sep{1.2, 1.2, 1.2};

            std::string name=calc_name({0,0,0}, dx, sep, dims);
            TestBase::add_test_factory(name, 0, [=](){
                return std::make_shared<TestMultipleMovingBeads>(name,  dx, sep, dims);
            });
        };

        vec3r_t v;
        for(v[0]=-1; v[0]<=1; v[0]++){
            for(v[1]=-1; v[1]<=1; v[1]++){
                for(v[2]=-1; v[2]<=1; v[2]++){
                    add(v, {7,8,9});
                    add(v, {4,8,8});
                    add(v, {4,8,4});
                }
            }
        }

        std::mt19937_64 rng(1);
        std::uniform_real_distribution<> urng(-1,1);
        for(int i=0; i<10; i++){
            add({
                urng(rng), urng(rng), urng(rng)
            }, {7,8,9});
            add({
                urng(rng), urng(rng), urng(rng)
            }, {4,4,4});
            add({
                urng(rng), urng(rng), urng(rng)
            }, {4,8,4});
            add({
                urng(rng), urng(rng), urng(rng)
            }, {4,8,8});
        }
    }

    TestMultipleMovingBeads(std::string name, const vec3r_t &v=vec3r_t{0.5,0,0}, const vec3r_t &sep=vec3r_t{1.2, 1.2, 1.2}, const vec3r_t &dims={7,8,9})
        : m_name(name)
        , m_v(v)
        , m_sep(sep)
        , m_dims(dims)
    {}

    virtual std::string name() const 
    {
        return m_name;
    }

    WorldState create_world() override
    {

        WorldStateBuilder b(m_dims);
        WorldState &s=b.data();

        s.t=0;
        s.dt=1.0/64;

        b.add_bead_type("W");
        b.add_polymer_type("W", {"W"}, {}, {});

        vec3r_t pos={0,0,0};
        while(pos[0] + m_sep[0] < m_dims[0]){
            while(pos[1] + m_sep[1] < m_dims[1]){
                while(pos[2] + m_sep[2] < m_dims[2]){
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
        for(unsigned i=0; i<m_x0.size(); i++){
            vec3r_t pos=s.beads.at(i).x;
            vec3r_t full_x=m_v * s.t * s.dt + m_x0[i];
            vec3r_t red_x=full_x.reduce_to_box(s.box);
            //std::cerr<<"Pos = "<<s.beads.at(0).x<<"\n";
            double err= (red_x - pos ).l2_norm();
            if(err > sqrt(m_steps_done) * 1e-4){
                std::stringstream acc;
                acc<<"Bead is in the wrong position at t="<<s.t*s.dt<<". Expected "<<red_x<<", got "<<pos<<", err="<<err<<", tol="<<sqrt(m_steps_done) * 1e-4;
                throw TestFailedException(acc.str());
            }
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
