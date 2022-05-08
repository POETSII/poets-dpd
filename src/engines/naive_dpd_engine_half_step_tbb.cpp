#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb.hpp"

bool naive_dpd_engine_half_step_tbb_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_half_step_tbb",
    [](){ return std::make_shared<NaiveDPDEngineHalfStepTBB>(); }
);