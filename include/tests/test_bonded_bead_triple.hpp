#ifndef test_bonded_bead_triple_hpp
#define test_bonded_bead_triple_hpp

#include "test_base.hpp"

#include "../dpd_state_builder.hpp"

#include <iostream>

class TestBondedBeadTriple
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;

    std::string m_name;
    vec3r_t m_x0; // Centre between beads
    vec3r_t m_dx; // Bead 0 starts at x0-2*dx, Bead 1 at x0, Bead 2 at x0+2*dx
    double m_kappa;
    double m_r0;

    double m_r_equilibrium;

    double m_sum_dist=0;
    double m_count_dist=0;

    int m_reversals=0;

public:
    static std::string calc_name(const vec3r_t &x0, const vec3r_t &dx, double kappa, double r0)
    {
        std::stringstream tmp;
        tmp<<"TestBondedBeadTriple[x0="<<x0<<";dx="<<dx<<";kappa="<<kappa<<";r0="<<r0<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [=](vec3r_t dx, double kappa, double r0)
        {
            std::string name=calc_name({4,4,4}, dx, kappa, r0);
            TestBase::add_test_factory(name, 2, [=](){
                return std::make_shared<TestBondedBeadTriple>(name, vec3r_t{4,4,4}, dx, kappa, r0);
            });
        };

        add({0.2,0,0}, 70, 0.6);

        add({0.1,0.1,0}, 30, 0.55);

        add({0.1,0.1,-0.1}, 100, 0.65);
    }

    TestBondedBeadTriple(const std::string &name, const vec3r_t &x0, const vec3r_t &dx, double kappa=1, double r0=0.5)
        : m_name(name)
        , m_x0(x0)
        , m_dx(dx)
        , m_kappa(kappa)
        , m_r0(r0)
    {
        require( (dx*2).l2_norm() < 0.99, "dx is too large, bond is broken." );
    }

    virtual std::string name() const 
    { return m_name; }

    WorldState create_world() override
    {
        WorldStateBuilder b( {8,8,8} );
        WorldState &s=b.data();
        s.t=0;
        s.dt=1.0/100;

        b.add_bead_type("A");
        b.add_bead_type("B");
        b.add_polymer_type("D", {"A","B", "A"}, { {m_kappa, m_r0, 0, 1}, {m_kappa, m_r0, 1, 2} }, {});

        double Cd=1;
        double Cc=1;
        b.set_interaction_strength("A", "B", Cc, Cd);
        b.set_interaction_strength("A", "A", Cc, Cd);
        b.set_interaction_strength("B", "B", Cc, Cd);

        m_r_equilibrium = (m_kappa*m_r0+Cc)/(m_kappa+Cc);
        require(m_r_equilibrium*2 > 1, "R equilibrium too small (cross-triple conservative interactions.");

        b.add_polymer("D", {{m_x0-m_dx*2}, {m_x0}, {m_x0+m_dx*2}});

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        require_close(s.dt * m_steps_done, s.t, "t is not dt*nSteps.");
        require_close(s.dt, 1/100.0, "dt is wrong.");

        double dist01=distance(s, s.beads[0].x, s.beads[1].x).l2_norm();
        require(dist01 < 0.99, "Bond has snapped.");

        double dist12=distance(s, s.beads[1].x, s.beads[2].x).l2_norm();
        require(dist12 < 0.99, "Bond has snapped.");
        
        /*
        std::cerr<<"t="<<s.t<<", dist="<<dist<<", equib="<<m_r_equilibrium<<", mean="<<m_sum_dist/m_count_dist<<"\n";
        std::cerr<<"  b0: x="<<s.beads[0].x<<", v="<<s.beads[0].v<<", f="<<s.beads[0].f<<"\n";
        std::cerr<<"  b0: x="<<s.beads[1].x<<", v="<<s.beads[1].v<<", f="<<s.beads[1].f<<"\n";
        */

        if(s.t > 10){
            m_sum_dist += dist12 + dist01;
            m_count_dist += 1;
        }

        if(s.t < 20){
            m_steps_done += m_step_dist;
            return m_step_dist;
        }else{
            //std::cerr<<"# t="<<s.t<<", dt="<<s.dt<<", dist="<<dist01+dist12<<", equib="<<2*m_r_equilibrium<<", mean="<<m_sum_dist/m_count_dist<<"\n";
            double mean_dist=m_sum_dist / m_count_dist;
            require( std::abs(mean_dist - 2*m_r_equilibrium) < 0.05, "Mean distance is wrong.");
            return 0;
        }
    }
};

#endif
