#ifndef test_base_hpp
#define test_base_hpp

#include "dpd/core/dpd_state.hpp"

#include "dpd/core/vec3.hpp"

#include <stdexcept>
#include <typeinfo> 
#include <functional>
#include <map>
#include <memory>
#include <tuple>
#include <random>

#include "dpd/core/dpd_state_validator.hpp"

#include "dpd/core/dpd_state_builder.hpp"

struct TestFailedException
    : std::runtime_error
{
    TestFailedException(const std::string &msg)
        : std::runtime_error(msg)
    {}
};

class TestBase
{
public:
    virtual ~TestBase()
    {}

    using test_instance_factory_t =std::function<std::shared_ptr<TestBase>()>;

    virtual std::string name() const =0;

    virtual WorldState create_world() =0;

    virtual unsigned get_advance_count(WorldState &s) =0;

    static void add_test_factory(const std::string &name, int complexity, test_instance_factory_t factory)
    {
        if(get_name_to_index().find(name)!=get_name_to_index().end()){
            throw std::runtime_error("Duplicate factory name.");
        }
        
        unsigned index=get_factories().size();

        get_factories().push_back({name,complexity,factory});
        get_name_to_index().insert({name,index});
    }

    using factory_tuple_t = std::tuple<std::string,int,test_instance_factory_t>;

    static std::vector<factory_tuple_t> get_test_factories_ordered_by_name()
    {
        std::vector<factory_tuple_t> tmp(get_factories());
        std::sort(tmp.begin(), tmp.end(), [](const factory_tuple_t &a, const factory_tuple_t &b){
            if(std::get<0>(a) < std::get<0>(b)) return true;
            if(std::get<0>(a) > std::get<0>(b)) return false;
            return (std::get<1>(a) < std::get<1>(b));
        });
        return tmp;
    }

    static std::vector<factory_tuple_t> get_test_factories_ordered_by_complexity()
    {
        std::vector<factory_tuple_t> tmp(get_factories());
        std::sort(tmp.begin(), tmp.end(), [](const factory_tuple_t &a, const factory_tuple_t &b){
            if(std::get<1>(a) < std::get<1>(b)) return true;
            if(std::get<1>(a) > std::get<1>(b)) return false;
            return (std::get<0>(a) < std::get<0>(b));
        });
        return tmp;
    }
protected:
    static std::map<std::string,unsigned> &get_name_to_index()
    {
        static std::map<std::string,unsigned> f;
        return f;
    }

    static std::vector<std::tuple<std::string,int,test_instance_factory_t>> &get_factories()
    {
        static std::vector<std::tuple<std::string,int,test_instance_factory_t>> f;
        return f;
    }

    bool contains(const WorldState &s, const vec3r_t &x)
    {
        bool res=true;
        for(unsigned i=0; i<3; i++){
            res &= (0 <= x[i] && x[i] < s.box[i]);
        }
        return res;
    }

    vec3r_t distance(const WorldState &s, const vec3r_t &from, const vec3r_t &to)
    {
        require(contains(s, from), "Point not in box.");
        require(contains(s, to), "Point not in box.");

        vec3r_t res;
        for(unsigned i=0; i<3; i++){
            double a=to[i]-from[i];
            double b=a+(s.box[i]);
            double c=a-(s.box[i]);
            if(std::abs(b) < std::abs(a)){
                a=b;
            }
            if(std::abs(c) < std::abs(a)){
                a=c;
            }
            res[i]=a;
        }
        return res;
    }

    void require(bool cond, const std::string &msg)
    {
        if(!cond){
            throw TestFailedException(msg);
        }
    }

    void require_close(double ref, double got, const char *msg)
    {
        if( std::abs(ref-got) > 1e-6 ){
            std::stringstream acc;
            acc<<msg<<", ref="<<ref<<", got="<<got;
            require(false, acc.str());
        }
    }

    void require_close(const vec3r_t &ref, const vec3r_t &got, const char *msg)
    {
        if( (ref-got).l2_norm() > 1e-6 ){
            std::stringstream acc;
            acc<<msg<<", ref="<<ref<<", got="<<got;
            require(false, acc.str());
        }
    }

    void require_close(const vec3r_t &ref, const vec3r_t &got, double tol, const char *msg)
    {
        if( (ref-got).l2_norm() > tol ){
            std::stringstream acc;
            acc<<msg<<", ref="<<ref<<", got="<<got<<", tol="<<tol;
            require(false, acc.str());
        }
    }

    static vec3r_t random_x(const vec3r_t &box, std::mt19937_64 &r)
    {
        vec3r_t res;
        std::uniform_real_distribution dist;
        for(unsigned i=0; i<3; i++){
            res[i] =  box[i] * dist(r);
        }
        return res;
    }

    static vec3r_t random_x(const WorldState &s, std::mt19937_64 &r)
    {
        return random_x(s.box, r);
    }
};

#endif