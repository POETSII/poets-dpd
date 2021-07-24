#ifndef test_bonded_bead_pair_hpp
#define test_bonded_bead_pair_hpp

#include "test_base.hpp"

#include "../dpd_state_builder.hpp"

#include <iostream>

class TestBondedBeadPair
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    std::string m_name;
    vec3r_t m_x0; // Centre between beads
    vec3r_t m_dx; // Bead 0 starts at x0+dx, Bead 1 at x0-dx
    double m_kappa;
    double m_r0;

    double m_r_equilibrium;

    double m_prev_dist;
    double m_prev_prev_dist;

    double m_sum_dist=0;
    double m_count_dist=0;

    int m_reversals=0;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx, double kappa, double r0)
    {
        std::stringstream tmp;
        tmp<<"TestBondedBeadPair[x0="<<x0<<";dx="<<dx<<";kappa="<<kappa<<";r0="<<r0<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [=](vec3r_t dx, double kappa, double r0)
        {
            std::string name=calc_name({2,2,2}, dx, kappa, r0);
            TestBase::add_test_factory(name, 2, [=](){
                return std::make_shared<TestBondedBeadPair>(name, vec3r_t{2,2,2}, dx, kappa, r0);
            });
        };

        add({0.2,0,0}, 10, 0.5);
        add({0.2,0,0}, 10, 0.4);
        add({0.2,0,0}, 10, 0.3);

        add({0.2,0.1,0}, 20, 0.5);
        add({0.2,0.1,0}, 20, 0.4);
        add({0.2,0.1,0}, 20, 0.3);

        add({-0.1,0.1,0.1}, 20, 0.5);
        add({-0.1,0.1,0.1}, 20, 0.4);
        add({-0.1,0.1,0.1}, 20, 0.3);
    }

    TestBondedBeadPair(const std::string &name, const vec3r_t &x0, const vec3r_t &dx, double kappa=1, double r0=0.5)
        : m_name(name)
        , m_x0(x0)
        , m_dx(dx)
        , m_kappa(kappa)
        , m_r0(r0)
    {
        require( dx.l2_norm() < 0.49, "dx is too large, bond is broken." );
    }

    virtual std::string name() const 
    { return m_name; }

    WorldState create_world() override
    {
        WorldStateBuilder b( {4,4,4} );
        WorldState &s=b.data();
        s.t=0;
        s.dt=1.0/10;

        b.add_bead_type("A");
        b.add_bead_type("B");
        b.add_polymer_type("D", {"A","B"}, { {m_kappa, m_r0, 0, 1} }, {});

        double Cd=1;
        double Cc=1;
        b.set_interaction_strength("A", "B", Cc, Cd);

        /* We have two forces:
            dx = xi - xj;
            r = |dx|
            Fc = (1-r) * C 
            Fb = kappa*(r-r0)

            So equilibrium at:
              kappa*(r-r0) = (1-r) * C
            Or:
              r=(kappa*r0+C)/(kappa+C)
        */

        m_r_equilibrium = (m_kappa*m_r0+Cc)/(m_kappa+Cc);

        b.add_polymer("D", {{m_x0+m_dx}, {m_x0-m_dx}});

        m_prev_dist=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();
        m_prev_prev_dist=m_prev_dist+0.1; // Create artificial minima

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        require_close(s.dt * m_steps_done, s.t, "t is not dt*nSteps.");

        require_close( s.beads[0].v, -s.beads[1].v, "Beads should be moving in opposite directions at same speed.");

        double dist=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();

        require(dist < 0.99, "Bond has snapped.");

        // Get distance from origin to bead
        vec3r_t tmp=s.beads[0].x - m_x0;
        if(tmp[0]<0 && m_dx[0]>0){ tmp=-tmp; }
        if(tmp[1]<0 && m_dx[1]>0){ tmp=-tmp; }
        if(tmp[2]<0 && m_dx[2]>0){ tmp=-tmp; }

        if(s.t > 0){
            require_close( normalise(m_dx), normalise(tmp), 1e-4*sqrt(s.t/s.dt), "Bead 0 should always be on the same line.");
        }
        
        /*
        std::cerr<<"t="<<s.t<<", dist="<<dist<<", equib="<<m_r_equilibrium<<", mean="<<m_sum_dist/m_count_dist<<"\n";
        std::cerr<<"  b0: x="<<s.beads[0].x<<", v="<<s.beads[0].v<<", f="<<s.beads[0].f<<"\n";
        std::cerr<<"  b0: x="<<s.beads[1].x<<", v="<<s.beads[1].v<<", f="<<s.beads[1].f<<"\n";
        */

        if(s.t > 10){
            m_sum_dist += dist;
            m_count_dist += 1;
        }

        m_prev_prev_dist=m_prev_dist;
        m_prev_dist=dist;

        if(s.t < 20){
            m_steps_done += 1;
            return 1;
        }else{
            //std::cerr<<"# t="<<s.t<<", dist="<<dist<<", equib="<<m_r_equilibrium<<", mean="<<m_sum_dist/m_count_dist<<"\n";
            double mean_dist=m_sum_dist / m_count_dist;
            require( std::abs(mean_dist - m_r_equilibrium) < 0.01, "Mean distance is wrong.");
            return 0;
        }
    }
};

#endif
