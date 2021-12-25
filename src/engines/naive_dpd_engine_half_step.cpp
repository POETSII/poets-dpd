#include "dpd/engines/naive/naive_dpd_engine_half_step.hpp"

bool naive_dpd_engine_half_step_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_half_step",
    [](){ return std::make_shared<NaiveDPDEngineHalfStep>(); }
);