#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb_v2.hpp"

bool naive_dpd_engine_half_step_tbb_v2_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_half_step_tbb_v2",
    [](){ return std::make_shared<NaiveDPDEngineHalfStepTBBV2>(); }
);