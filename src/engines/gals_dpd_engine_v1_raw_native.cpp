#include "dpd/engines/gals/gals_dpd_engine_v1_raw.hpp"

bool gals_dpd_engine_v1_raw_native_registered = DPDEngineFactory::RegisterFactory(
    "gals_dpd_engine_v1_raw_native",
    [](){ return std::make_shared<GALSDPDEngineV1Raw>(); }
);