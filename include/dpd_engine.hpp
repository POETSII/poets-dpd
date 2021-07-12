#ifndef dpd_engine_hpp
#define dpd_engine_hpp

#include "dpd_state.hpp"

class DPDEngine
{
public:
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
};

#endif