#ifndef test_isolated_bead_hpp
#define test_isolated_bead_hpp

#include "test_base.hpp"


#include <iostream>

class TestCollidingBeadPair
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

    int m_reversals=0;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx)
    {
        std::stringstream tmp;
        tmp<<"TestCollidingBeadPair[x0="<<x0<<";dx="<<dx<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [](vec3r_t dx)
        {
            vec3r_t x0{2,2,2};
            std::string name=calc_name(x0, dx);
            TestBase::add_test_factory(name, 1, [=](){
                return std::make_shared<TestCollidingBeadPair>(name, x0, dx);
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

    }

    TestCollidingBeadPair(const std::string &name, const vec3r_t &x0, const vec3r_t &dx)
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
        s.interactions[0].dissipative=0; // Unit conservative, zero dissipative

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
        require_close(m_steps_done, s.t, "t is not dt*nSteps.");

        double dist=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();

        if(m_dx.l2_norm() > 0.5){
            //std::cerr<<"dx="<<m_dx<<"\n;";
            require_close( m_x0+m_dx,  s.beads[0].x, "Bead 0 should not move.");
            require_close( m_x0-m_dx,  s.beads[1].x, "Bead 1 should not move.");
        }else{
            // Get distance from origin to bead
            vec3r_t tmp=s.beads[0].x - m_x0;
            if(tmp[0]<0 && m_dx[0]>0){ tmp=-tmp; }
            if(tmp[1]<0 && m_dx[1]>0){ tmp=-tmp; }
            if(tmp[2]<0 && m_dx[2]>0){ tmp=-tmp; }
            // Low tolerance, as they drift quite quickly in single precision
            if(s.t>0){
                require_close( normalise(m_dx), normalise(tmp), 1e-3*sqrt(s.t/s.dt), "Bead 0 should always be on the same line.");
            }
            
            if( s.t == 2){
                require_close( normalise(m_dx), normalise(s.beads[0].v), 1e-5, "Bead 0 should move along +dx at time-step 1." );
                require_close( normalise(-m_dx), normalise(s.beads[1].v), 1e-5, "Bead 1 should move along -dx at time-step 1." );
            }

            if(m_prev_dist < 1 && dist < 1){
                if(m_prev_dist < dist){
                    require( m_prev_v[0].l2_norm() < s.beads[0].v.l2_norm(), "Bead 0 should be accelerating.");
                }else if(m_prev_dist > dist){
                     require( m_prev_v[0].l2_norm() > s.beads[0].v.l2_norm(), "Bead 0 should be decelerating.");
                }
            }

            if(m_prev_dist < dist && m_prev_dist < m_prev_prev_dist ) {
                m_reversals++;
                // there is no energy loss
                //std::cerr<<"# Reversal at t="<<s.t-s.dt<<", x[0]="<<m_prev_x[0]<<", x[1]"<<m_prev_x[1]<<", dist="<<m_prev_dist<<", distP="<<m_prev_prev_dist<<", distN="<<dist<<"\n";
                require( std::abs(m_dx.l2_norm()*2 - m_prev_dist) < 0.001, "Closest distance should be very close to dx." );
            }

            if(dist < 1.1){
                //std::cerr<<"t="<<s.t<<", x0="<<s.beads[0].x<<", x1="<<s.beads[1].x<<", dist="<<dist<<", v="<<s.beads[0].v<<", f="<<s.beads[0].f<<"\n";
            }
            require_close( s.beads[0].v, -s.beads[1].v, "Beads should be moving in opposite directions at same speed.");
        }

        m_prev_v[0]=s.beads[0].v;
        m_prev_v[1]=s.beads[1].v;
        m_prev_prev_dist=m_prev_dist;
        m_prev_dist=dist;

        m_prev_x[0]=s.beads[0].x;
        m_prev_x[1]=s.beads[1].x;

        if(m_reversals < 4){
            m_steps_done += m_step_dist;
            return m_step_dist;
        }else{
            return 0;
        }
    }
};

#endif
