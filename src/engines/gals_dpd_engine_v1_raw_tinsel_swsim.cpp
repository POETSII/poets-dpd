#include "dpd/engines/gals/gals_dpd_engine_v1_raw_tinsel.hpp"

bool gals_dpd_engine_v1_raw_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "gals_dpd_engine_v1_raw_tinsel_swsim",
    [](){ return std::make_shared<GALSDPDEngineV1RawTinsel<POLiteSWSim<>>>(); }
);