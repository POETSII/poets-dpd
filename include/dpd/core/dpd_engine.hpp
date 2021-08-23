#ifndef dpd_engine_hpp
#define dpd_engine_hpp

#include "dpd_state.hpp"

#include <memory>
#include <map>
#include <functional>

class DPDEngine
{
public:
    /* Check the given world and check it can be attached to this engine.
        Return empty string if ok, otherwise return the reason it can't
        be attached.
    */
    virtual std::string CanSupport(const WorldState *) const
    { return ""; }

    // Attach the given world-state to this engine.
    // Any future methods are relative to this state.
    // Attach(nullptr) will detach the engine.
    // While attached, the engine can assume that all pointers
    // are stable (e.g. no vectors moved), the types will all be
    // constant, and the bead states will not be changed.
    virtual void Attach(WorldState *state) =0;

    // Run the worldstate forwards for the given number
    // of steps. While executing the engine has exclusive
    // read-write access to the beads array mutable properties,
    // but it must not re-order or otherwise change the vector.
    // An engine can reasonably assume that nSteps will be large
    // enough to ammortise any setup costs.
    virtual void Run(unsigned nSteps) =0;

    // Run the world state for a number of intervals of given
    // size, with the callback called at the end of each interval.
    // When the callback is invoked the world-state will be valid,
    // though execution may still be proceeding in parallel.
    // The callback should return true to continue execution, or
    // false to terminate early
    // Returns the number of steps actually completed.
    // If interval_count is -1 then run for ever.
    virtual unsigned Run(
        int interval_count,
        unsigned interval_size,
        std::function<bool()> interval_callback
    ) {
        unsigned done=0;
        while(1){
            if(interval_count==0){
                break;
            }
            Run(interval_size);
            done += interval_size;
            if(!interval_callback()){
                break;
            }
            if(interval_count>0){
                --interval_count;
            }
        }
        return done;
    }
};

class DPDEngineFactory
{
public:
    using creator_func_t = std::function<std::shared_ptr<DPDEngine>()>;

    static std::vector<std::string> ListFactories()
    {
        auto &m=get_map();
        std::vector<std::string> res;
        for(const auto &e : m){
            res.push_back(e.first);
        }
        return res;
    }
    
    static bool RegisterFactory(const std::string &name, creator_func_t func)
    {
        auto &m=get_map();
        auto it=m.find(name);
        if(it!=m.end()){
            throw std::runtime_error("Duplicate DPD engines registered as '"+name+"'");
        }
        m.insert({name,func});
        return true;
    }

    static std::shared_ptr<DPDEngine> create(const std::string &name)
    {
        auto &m=get_map();
        auto it=m.find(name);
        if(it==m.end()){
            throw std::runtime_error("No DPD engine registered called '"+name+"'");
        }
        return it->second();
    }
private:
    using factory_map=std::map<std::string,creator_func_t>;

    static factory_map &get_map()
    {
        static factory_map the_map;
        return the_map;
    }
};

#endif