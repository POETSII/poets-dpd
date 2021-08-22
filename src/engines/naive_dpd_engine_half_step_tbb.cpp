#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb.hpp"

bool naive_dpd_engine_half_step_tbb_registered = DPDEngineFactory::RegisterFactory(
    "NaiveDPDEngineHalfStepTBB",
    [](){ return std::make_shared<NaiveDPDEngineHalfStepTBB>(); }
);