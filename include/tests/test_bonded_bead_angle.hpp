#ifndef test_bonded_bead_angle_hpp
#define test_bonded_bead_angle_hpp

#include "test_base.hpp"

#include "../dpd_state_builder.hpp"

#include <iostream>

class TestBondedBeadAngle
    : public TestBase
{
    unsigned m_step_dist=1;
    unsigned m_steps_done=0;
    

    std::string m_name;
    vec3r_t m_x1; // Centre between beads
    vec3r_t m_dx0; // Bead 0 starts at x1+dx0
    vec3r_t m_dx2; // Bead 2 starts at x1+dx1
    double m_kappa;
    double m_r0;
    double m_theta0;

    double m_r_equilibrium;

    double m_sum_dist=0;
    double m_sum_angle=0;
    double m_count=0;

public:
    static std::string calc_name(const vec3r_t &x1, const vec3r_t &dx0, const vec3r_t &dx2, double kappa, double r0, double theta0)
    {
        std::stringstream tmp;
        tmp<<"TestBondedBeadAngle[x1="<<x1<<";dx0="<<dx0<<";dx2="<<dx2<<";kappa="<<kappa<<";r0="<<r0<<";theta0="<<theta0<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [=](vec3r_t x0, vec3r_t dx0, vec3r_t dx2, double kappa, double r0, double theta0)
        {
            std::string name=calc_name(x0, dx0, dx2, kappa, r0, theta0);
            TestBase::add_test_factory(name, 3, [=](){
                return std::make_shared<TestBondedBeadAngle>(name, x0, dx0, dx2, kappa, r0, theta0);
            });
        };

        std::mt19937_64 rng;

        static const vec3r_t box{8,8,8};

        auto rp=[&]() -> vec3r_t { return random_x(box, rng); };

        add({0,0,0}, {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.0);
        add({0,0,0}, {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.1);
        add({0,0,0}, {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.2);

        add({0,0,0}, {-0.4,0,0}, {0.3,0.2,0}, 40, 0.50, 0.0);
        add({0,0,0}, {-0.4,0,0}, {0.3,0.2,0}, 40, 0.50, 0.1);
        add({0,0,0}, {-0.4,0,0}, {0.3,0.1,0}, 40, 0.50, 0.2);

        add({0,0,0}, {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.0);
        add({0,0,0}, {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.1);
        add({0,0,0}, {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.2);


        add(rp(), {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.0);
        add(rp(), {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.1);
        add(rp(), {-0.4,0,0}, {0.4,0.1,0}, 40, 0.50, 0.2);

        add(rp(), {-0.4,0,0}, {0.3,0.2,0}, 40, 0.50, 0.0);
        add(rp(), {-0.4,0,0}, {0.3,0.2,0}, 40, 0.50, 0.1);
        add(rp(), {-0.4,0,0}, {0.3,0.1,0}, 40, 0.50, 0.2);

        add(rp(), {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.0);
        add(rp(), {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.1);
        add(rp(), {-0.4,0,0}, {0.2,0.1,0.2}, 40, 0.50, 0.2);


    }

    TestBondedBeadAngle(const std::string &name, const vec3r_t &x1, const vec3r_t &dx0, const vec3r_t &dx2, double kappa=1, double r0=0.5, double theta0=0)
        : m_name(name)
        , m_x1(x1)
        , m_dx0(dx0)
        , m_dx2(dx2)
        , m_kappa(kappa)
        , m_r0(r0)
        , m_theta0(theta0)
    {
        require( (dx0).l2_norm() < 0.99, "dx0 is too large, bond is broken." );
        require( (dx2).l2_norm() < 0.99, "dx2 is too large, bond is broken." );
    }

    virtual std::string name() const 
    { return m_name; }

    WorldState create_world() override
    {
        static const vec3r_t box{8,8,8};

        WorldStateBuilder b( box);
        WorldState &s=b.data();
        s.t=0;
        s.dt=0.005;

        b.add_bead_type("A");
        b.add_bead_type("B");
        b.add_polymer_type("D",
            {"A","B", "A"},
            { {m_kappa, m_r0, 0, 1}, {m_kappa, m_r0, 1, 2} },
            { {m_kappa, m_theta0, 0, 1} }
        );

        double Cd=1;
        double Cc=1;
        b.set_interaction_strength("A", "B", Cc, Cd);
        b.set_interaction_strength("A", "A", Cc, Cd);
        b.set_interaction_strength("B", "B", Cc, Cd);

        m_r_equilibrium = (m_kappa*m_r0+Cc)/(m_kappa+Cc);
        require(m_r_equilibrium*2 > 1, "R equilibrium too small (cross-triple conservative interactions.");

        b.add_polymer("D", {{m_x1+m_dx0}, {m_x1}, {m_x1+m_dx2}}, true);

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        require_close(s.dt * m_steps_done, s.t, "t is not dt*nSteps.");

        auto dx01=distance(s, s.beads[0].x, s.beads[1].x);
        auto dx12=distance(s, s.beads[1].x, s.beads[2].x);

        double dist01=dx01.l2_norm();
        require(dist01 < 0.99, "Bond has snapped.");

        double dist12=dx12.l2_norm();
        require(dist12 < 0.99, "Bond has snapped.");
        
        
        //std::cout<<s.t<<", "<<angle(dx01,dx12)<<", "<<m_theta0<<"\n";
    
        //std::cerr<<"t="<<s.t<<", dist1="<<dist01<<", dist12="<<dist12<<", equib="<<m_r_equilibrium<<", mean="<<m_sum_dist/m_count<<"\n";
        //std::cerr<<"  b0: x="<<s.beads[0].x<<", v="<<s.beads[0].v<<", f="<<s.beads[0].f<<"\n";
        //std::cerr<<"  b1: x="<<s.beads[1].x<<", v="<<s.beads[1].v<<", f="<<s.beads[1].f<<"\n";
        //std::cerr<<"  b2: x="<<s.beads[2].x<<", v="<<s.beads[2].v<<", f="<<s.beads[2].f<<"\n";
        

        if(s.t > 10){
            m_sum_dist += dist12 + dist01;
            double a=angle(dx01, dx12);
            m_sum_angle += a;
            require(std::isnormal(m_sum_angle), "sum_angle is not normal.");
            m_count += 1;
        }

        // Take some energy out of the system, as there are no
        // other beads to absorb
        for(unsigned i=0; i<3; i++){
            s.beads[i].v = s.beads[i].v * 0.99;
        }

        if(s.t < 100){
            m_steps_done += 1;
            return 1;
        }else{
            //std::cerr<<"# t="<<s.t<<", dist="<<dist01+dist12<<", equib="<<2*m_r_equilibrium<<", mean="<<m_sum_dist/m_count<<"\n";
            //std::cerr<<"# t="<<s.t<<", angle="<<angle(dx01,dx12)<<", target="<<m_theta0<<", mean="<<m_sum_angle/m_count<<"\n";
            double mean_dist=m_sum_dist / m_count;
            require( std::abs(mean_dist - 2*m_r_equilibrium) < 0.1, "Mean distance is wrong.");
            double mean_angle=m_sum_angle / m_count;
            require( std::abs(mean_angle - m_theta0) < 0.1, "Mean angle is wrong.");
            return 0;
        }
    }
};

#endif
