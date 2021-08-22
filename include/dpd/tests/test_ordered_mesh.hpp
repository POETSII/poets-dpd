#ifndef test_ordered_mesh_hpp
#define test_ordered_mesh_hpp

#include "test_base.hpp"

#include "dpd/core/dpd_state_to_vtk.hpp"

#include <iostream>

class TestOrderedMesh
    : public TestBase
{
    std::string m_name;
    vec3i_t m_dims;
public:
    static std::string calc_name(const vec3i_t &dims)
    {
        std::stringstream tmp;
        tmp<<"TestOrderedMesh[dims="<<dims<<"]";
        return tmp.str();
    }

    static void register_tests()
    {
        auto add = [=](vec3i_t dims)
        {
            std::string name=calc_name(dims);
            TestBase::add_test_factory(name, 3, [=](){
                return std::make_shared<TestOrderedMesh>(name, dims);
            });
        };

        add({4,4,4});
    }

    TestOrderedMesh(const std::string &name, const vec3i_t &dims)
        : m_name(name)
        , m_dims(dims)
    {
    }

    virtual std::string name() const 
    { return m_name; }

    WorldState create_world() override
    {
        WorldStateBuilder b{vec3r_t(m_dims)};
        WorldState &s=b.data();
        s.t=0;
        s.dt=1.0/100;

        b.add_bead_type("A");
        b.add_polymer_type("A", {"A"}, {}, {});

        b.set_interaction_strength(0, 0, 25, 4.5);

        int beads_per_dim=6;

        double density=1;
        for(unsigned i=0; i<3; i++){
            density *= beads_per_dim;
            density /= m_dims[i];
        }
        assert(density <= 4);

        vec3i_t pos;
        for(pos[0]=0; pos[0]<beads_per_dim; pos[0]++){
            for(pos[1]=0; pos[1]<beads_per_dim; pos[1]++){
                for(pos[2]=0; pos[2]<beads_per_dim; pos[2]++){
                    vec3r_t x;
                    for(int i=0; i<3; i++){
                        x[i] = pos[i] / double(beads_per_dim) * m_dims[i];
                    }
                    b.add_polymer("A", {{x}});
                    //std::cerr<<"x="<<x<<"\n";
                }
            }
        }

        return b.extract();
    }

    unsigned get_advance_count(WorldState &s) override
    {
        double temperature=0;
        vec3r_t sum_v;
        vec3r_t sum_f;
        for(const auto &b : s.beads){
            temperature += b.v.dot(b.v);
            sum_v += b.v;
            sum_f += b.f;
        }
        temperature /= (s.beads.size()*3);

        double density=s.beads.size() / double(m_dims[0]*m_dims[1]*m_dims[2]);

        if(0 && std::abs(round(s.t)-s.t) < 1e-9 ){
            std::cerr<<"t="<<s.t<<", KbT="<<temperature<<", density="<<density<<"\n";
            std::cerr<<"  v[0]="<<s.beads[0].v<<", mean_v="<<sum_v*(1.0/s.beads.size())<<", mean_f="<<sum_f*(1.0/s.beads.size())<<"\n";
        }

        // TODO : What should the temperature actually be. It doesnt seem to be 1...
        
        if(s.t < 10){
            return 1;
        }else{
            return 0;
        }
    }
};

#endif
