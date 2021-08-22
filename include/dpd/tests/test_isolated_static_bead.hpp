#ifndef test_isolated_bead_hpp
#define test_isolated_bead_hpp

#include "test_base.hpp"

class TestIsolatedStaticBead
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    virtual std::string name() const 
    { return "IsolatedStaticBead"; }

    WorldState create_world() override
    {
        WorldState res;
        for(unsigned i=0; i<3; i++){
            res.origin[i]=0;
            res.box[i]=2;
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
            { 0, 0, 0 },
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
        if(pos != vec3r_t()){
            throw TestFailedException("Bead moved away from zero.");
        }
        if( abs(s.t - s.dt * m_steps_done ) > 1e-10 ){
            throw TestFailedException("t is not dt*nSteps.");
        }
        if(s.t < 1){
            m_steps_done += m_step_dist;
            return m_step_dist++;
        }else{
            return 0;
        }
    }
};

#endif
