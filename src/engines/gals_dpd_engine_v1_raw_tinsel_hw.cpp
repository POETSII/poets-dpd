#include "dpd/engines/gals/gals_dpd_engine_v1_raw_tinsel.hpp"

#include "POLiteHW.h"

bool gals_dpd_engine_v1_raw_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "gals_dpd_engine_v1_raw_tinsel_hw",
    [](){ return std::make_shared<GALSDPDEngineV1RawTinsel<POLiteHW<>>>(); }
);