#ifndef dpd_engine_hpp
#define dpd_engine_hpp

#include "dpd_state.hpp"

#include <memory>
#include <map>
#include <functional>

/*
class EngineLog
{
public:
    virtual ~EngineLog()
    {}

    virtual bool is_active() const
    { return false; }

    virtual void log_end_step(double t)
    {}

    virtual void log_pairwise_force(uint32_t hash1, uint32_t hash2, const vec3r_t &dx, double r, const vec3r_t &f, double fCon, double fDiss, double fRan, double fBond)
    {}

    virtual void log_angle_force(uint32_t hashHead, uint32_t hashMid, uint32_t hashTail, const vec3r_t &fHead, const vec3r_t &fMid, const vec3r_t &fTail)
    {}

    virtual void log_bead(uint32_t hash, const vec3r_t &x, const vec3r_t &v, const vec3_t &f)
    {}

    virtual void log_end_step(double t)
    {}
};


class StreamEngineLog
    : public EngineLog
{
private:
    std::ostream &dst;
    double m_t;
public:
    virtual bool is_active() const
    { return true; }

    void log_begin_step(double t)
    {
        m_t=t;
    }

    void log_pairwise_force(uint32_t hash1, uint32_t hash2, const vec3r_t &dx, double r, const vec3r_t &f, double fCon, double fDiss, double fRan, double fBond)
    {
        dst << "Pair,"<<m_t<<","<<hash1<<","<<hash2<<","<<dx<<","<<r<<","<<f<<","<<fCon<<","<<fDiss<<","<<fRan<<","<<fBond<<"\n";
    }

    void log_angle_force(uint32_t hashHead, uint32_t hashMid, uint32_t hashTail, const vec3r_t &fHead, const vec3r_t &fMid, const vec3r_t &fTail)
    {
        dst << "Angle,"<<m_t<<","<<hashHead<<","<<hashMid<<","<<hashTail<<","<<fHead<<","<<fMid<<","<<fTail<<"\n";
    }

    void log_bead_end_step(uint32_t hash, const vec3r_t &x, const vec3r_t &v, const vec3r_t &f)
    {
        dst << "BeadEnd,"<<m_t<<","<<hash<<","<<x<<","<<v<<","<<f<<"\n";
    }

    void log_end_step(double t)
    {
        assert(t==m_t);
    }
};
*/

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