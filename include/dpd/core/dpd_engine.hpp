#ifndef dpd_engine_hpp
#define dpd_engine_hpp

#include "dpd_state.hpp"

#include <memory>
#include <map>
#include <functional>

class DPDEngine
{
protected:
    virtual bool CanSupportHookeanBonds() const
    { return false; }

    virtual bool CanSupportAngleBonds() const
    { return false; }

    virtual bool CanSupportStationaryBeadTypes() const
    { return false; }

    virtual bool CanSupportSpatialDPDParameters() const
    { return false; }

    virtual bool CanSupportGeneralDPDParameters() const
    { return false; }

    virtual bool CanSupportSpatialBondParameters() const
    { return false; }

    virtual bool CanSupportGeneralBondParameters() const
    { return false; }

public:
    virtual ~DPDEngine()
    {}

    struct timings_t
    {
        double compile = -1; // Everything todo with translating the in-memory WorldState into a loadable graph
        double aquire = -1;  // Any extra time taken to aquire resources
        double configure = -1; // Time taken to move loadable graph into hardware
        double execute_to_first_bead = -1; // Time taken from start of execution till first bead in first batch
        double execute_to_last_bead = -1;  // Time taken from start of execution till last bead in last batch
        double perf_counters = -1;         // Any extra time needed to extract performance counters
    };


    /* Check the given world and check it can be attached to this engine.
        Return empty string if ok, otherwise return the reason it can't
        be attached.

        Engines are recommended to call this as the base and override
        like CanSupportXXX(), so that they can reject things they didn't
        know about at time of writing.
    */
    virtual std::string CanSupport(const WorldState *s) const
    {
        std::unordered_set<std::string> space_params{"x","y","z", "ux", "uy", "uz"};

        std::unordered_set<std::string> interaction_vars;
        for(const auto &ii : s->interactions){
            ii.conservative.collect_variables(interaction_vars);
            ii.dissipative.collect_variables(interaction_vars);
        }
        if(!interaction_vars.empty()){
            bool space_only=true;
            for(auto k : interaction_vars){
                if(space_params.find(k)==space_params.end()){
                    space_only=false;
                    break;
                }
            }
            if(space_only){
                if(!CanSupportSpatialDPDParameters()){
                    return "Engine does not support space varying DPD parameters.";
                }
            }else{
                if(!CanSupportGeneralDPDParameters()){
                    return "Engine does not support aribtrary varying DPD parameters.";
                }
            }
        }

        for(auto &b : s->bead_types){
            if(b.stationary && !CanSupportStationaryBeadTypes()){
                return "Engine does not support stationary bead types.";
            }
        }

        std::unordered_set<std::string> bond_vars;
        for(const auto &p : s->polymer_types){
            if(!p.bonds.empty() && !CanSupportHookeanBonds()){
                return "Engine does not support hookean bonds.";
            }
            if(!p.bond_pairs.empty() && !CanSupportAngleBonds()){
                return "Engine does not support angle bonds.";
            }

            for(auto b : p.bonds){
                b.kappa.collect_variables(bond_vars);
                b.r0.collect_variables(bond_vars);
            }
            for(auto bp : p.bond_pairs){
                bp.kappa.collect_variables(bond_vars);
                bp.theta0.collect_variables(bond_vars);
            }
        }
        if(!bond_vars.empty()){
            bool space_only=true;
            for(auto k : bond_vars){
                if(space_params.find(k)==space_params.end()){
                    space_only=false;
                    break;
                }
            }
            if(space_only){
                if(!CanSupportSpatialBondParameters()){
                    return "Engine does not support space varying bond parameters.";
                }
            }else{
                if(!CanSupportGeneralBondParameters()){
                    return "Engine does not support aribtrary varying bond parameters.";
                }
            }
        }
        return "";
    }

    virtual double GetMaxBondLength() const
    { return 1.0; }

    /* If the engine needs to aquire hardware, this tells it to start
        aquiring resources in the background (if it can).
    */
    virtual void PrepareResources()
    {}

    /*! If the system supports timings, this will export them and return true. Otherwise it returns false.
        Timings are really only intended to support bench-mark style executions, so
        something like:
         Attach();
         Run();
    */
    virtual bool GetTimings(timings_t &timings)
    {
        return false;
    }

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