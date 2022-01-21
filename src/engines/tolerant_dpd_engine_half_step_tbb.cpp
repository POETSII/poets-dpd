#include "dpd/engines/naive/tolerant_dpd_engine_half_step_tbb.hpp"

bool tolerant_dpd_engine_half_step_tbb_registered = DPDEngineFactory::RegisterFactory(
    "tolerant_dpd_engine_half_step_tbb",
    [](){ return std::make_shared<TolerantDPDEngineHalfStepTBB>(); }
);