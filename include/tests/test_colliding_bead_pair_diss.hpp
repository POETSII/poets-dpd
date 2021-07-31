#ifndef test_colliding_bead_pair_diss_hpp
#define test_colliding_bead_pair_diss_hpp

#include "test_base.hpp"

#include "../dpd_state_builder.hpp"

#include <iostream>
#include <climits>
#include <cfloat>

class TestCollidingBeadPairDiss
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    std::string m_name;
    vec3r_t m_x0; // Centre between beads
    vec3r_t m_dx; // Bead 0 at x0+dx, Bead 1 at x0-dx

    vec3r_t m_prev_v[2];
    vec3r_t m_prev_x[2];
    double m_prev_dist;
    double m_prev_prev_dist;

    double m_last_closest=0;

    int m_reversals=0;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx)
    {
        std::stringstream tmp;
        tmp<<"TestCollidingBeadPairDiss[x0="<<x0<<";dx="<<dx<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [](vec3r_t dx, vec3r_t x0={2,2,2})
        {
            std::string name=calc_name(x0, dx);
            TestBase::add_test_factory(name, 1, [=](){
                return std::make_shared<TestCollidingBeadPairDiss>(name, x0, dx);
            });
        };

        add({0.2,0,0});
        add({0,0.2,0});
        add({0,0,0.2});

        add({-0.2,0,0});
        add({0,-0.2,0});
        add({0,0,-0.2});

        add({0.2,0.2,0});
        add({0,0.2,0.2});
        add({0.2,0,0.2});

        add({-0.2,0.2,0});
        add({0,-0.2,0.2});
        add({0.2,0,-0.2});
        

        add({0.15,0.15,0.15});
        add({0.15,-0.15,0.15});
        add({0.15,-0.15,-0.15});
        add({-0.15,-0.15,-0.15});

        add({0.2,0,0}, {1,1,1});
        add({0,0.2,0}, {1,2,3});
        add({0,0,0.2}, {1.5,1.5,1.5});
    }

    TestCollidingBeadPairDiss(const std::string &name, const vec3r_t &x0, const vec3r_t &dx)
        : m_name(name)
        , m_x0(x0)
        , m_dx(dx)
    {}

    virtual std::string name() const 
    { return m_name; }

    WorldState create_world() override
    {
        WorldStateBuilder b( {4,4,4} );
        WorldState &s=b.data();
        s.t=0;
        s.dt=1.0/256;

        b.add_bead_type("W");
        b.add_polymer_type("W", {"W"}, {}, {});

        s.interactions[0].conservative=1;
        s.interactions[0].dissipative=0.1; // Small dissipative

        b.add_monomer("W", m_x0+m_dx);
        b.add_monomer("W", m_x0-m_dx);

        m_prev_v[0]=vec3r_t();
        m_prev_v[1]=vec3r_t();
        m_prev_dist=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();
        m_prev_prev_dist=m_prev_dist+0.1; // Create artificial minima

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        require_close(s.dt * m_steps_done, s.t, "t is not dt*nSteps.");

        double dist=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();

        if(m_dx.l2_norm() > 0.5){
            //std::cerr<<"dx="<<m_dx<<"\n;";
            require_close( m_x0+m_dx,  s.beads[0].x, "Bead 0 should not move.");
            require_close( m_x0-m_dx,  s.beads[1].x, "Bead 1 should not move.");
        }else if(s.t>s.dt){
            // Get distance from origin to bead
            vec3r_t tmp=distance(s, m_x0, s.beads[0].x);
            // Low tolerance, as they drift quite quickly in single precision
            //std::cerr<<"  line="<<tmp<<", norm="<<normalise(tmp)<<", ref="<<normalise(m_dx)<<", diff="<< (normalise(m_dx)-normalise(tmp)).l2_norm() <<"\n";
            require_close( normalise(m_dx), normalise(tmp), 1e-3*sqrt(s.t/s.dt), "Bead 0 should always be on the same line.");

            if( std::abs(s.t-s.dt) < 1e-9){
                require_close( normalise(m_dx), normalise(s.beads[0].v), "Bead 0 should move along +dx at time-step 1." );
                require_close( normalise(-m_dx), normalise(s.beads[1].v), "Bead 1 should move along -dx at time-step 1." );
            }else if(std::abs(s.t-2*s.dt) < 1e-9){
            }

            if(m_prev_dist < dist && m_prev_dist < m_prev_prev_dist ) {
                m_reversals++;
                //std::cerr<<"# Reversal at t="<<s.t-s.dt<<", x[0]="<<m_prev_x[0]<<", x[1]"<<m_prev_x[1]<<", dist="<<m_prev_dist<<", distP="<<m_prev_prev_dist<<", distN="<<dist<<"\n";
                m_last_closest=m_prev_dist;
            }

            if(dist < 1.1){
                //std::cerr<<"t="<<s.t<<", x0="<<s.beads[0].x<<", x1="<<s.beads[1].x<<", dist="<<dist<<", v="<<s.beads[0].v<<", f="<<s.beads[0].f<<"\n";
            }
            require_close( s.beads[0].v, -s.beads[1].v, "Beads should be moving in opposite directions at same speed.");
        
            require( s.beads[0].v.l2_norm()>1e-6, "Beads should not stop moving.");
        }

        m_prev_v[0]=s.beads[0].v;
        m_prev_v[1]=s.beads[1].v;
        m_prev_prev_dist=m_prev_dist;
        m_prev_dist=dist;

        m_prev_x[0]=s.beads[0].x;
        m_prev_x[1]=s.beads[1].x;

        /* Note: this test is really sensitive to rounding errors in the bead position.
            Once they move slightly off-axis, they tend to push each other further off,
            so in single-precision it will be good for the early stages, then rapidly
            move off the original dx axis.

            Note that if an implementation stores absolute x coordinates as floating-point then 
            asymmetries will build up due to quantisation, so eventually all implementations
            will see this behaviour.
        */
        if(m_reversals < 3){
            m_steps_done += m_step_dist;
            return 1;
        }else{
            return 0;
        }
    }
};

#endif
