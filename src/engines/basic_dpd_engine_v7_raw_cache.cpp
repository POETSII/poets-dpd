#include "dpd/engines/basic/basic_dpd_engine_v7_raw.hpp"

bool basic_dpd_engine_v7_raw_cache_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v7_raw_cache",
    [](){ return std::make_shared<BasicDPDEngineV7Raw<true>>(); }
);