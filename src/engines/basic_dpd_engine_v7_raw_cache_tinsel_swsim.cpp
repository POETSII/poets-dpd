#include "dpd/engines/basic/basic_dpd_engine_v7_raw_tinsel.hpp"

bool basic_dpd_engine_v7_raw_cache_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v7_raw_cache_tinsel_swsim",
    [](){ return std::make_shared<BasicDPDEngineV7RawTinsel<POLiteSWSim<>,true>>(); }
);