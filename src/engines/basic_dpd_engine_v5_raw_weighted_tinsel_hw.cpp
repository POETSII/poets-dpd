#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

bool basic_dpd_engine_v5_raw_tinsel_weighted_hw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v5_raw_weighted_tinsel_hw",
    [](){ return std::make_shared<BasicDPDEngineV5RawTinsel<POLiteHW<>>>(true); }
);