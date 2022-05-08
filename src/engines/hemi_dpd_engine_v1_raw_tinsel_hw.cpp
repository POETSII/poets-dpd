#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw_tinsel.hpp"

#include "POLiteHW.h"

bool hemi_dpd_engine_v1_raw_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "hemi_dpd_engine_v1_raw_tinsel_hw",
    [](){ return std::make_shared<HemiDPDEngineV1RawTinsel<POLiteHW<>>>(); }
);