#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw_tinsel.hpp"

#include "POLiteSWSim_PGraph.h"

bool hemi_dpd_engine_v1_raw_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "hemi_dpd_engine_v1_raw_tinsel_swsim",
    [](){ return std::make_shared<HemiDPDEngineV1RawTinsel<POLiteSWSim<>>>(); }
);