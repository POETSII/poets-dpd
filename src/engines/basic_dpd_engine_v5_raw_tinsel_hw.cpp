#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

bool basic_dpd_engine_v5_raw_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV5RawTinselHW",
    [](){ return std::make_shared<BasicDPDEngineV5RawTinsel<POLiteHW<>>>(); }
);

bool basic_dpd_engine_v5_raw_tinsel_weighted_hw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV5RawTinselWeightedHW",
    [](){ return std::make_shared<BasicDPDEngineV5RawTinsel<POLiteHW<>>>(true); }
);