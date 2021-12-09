#include "dpd/engines/basic/basic_dpd_engine_v7_raw.hpp"

bool basic_dpd_engine_v7_raw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV7Raw",
    [](){ return std::make_shared<BasicDPDEngineV7Raw<false>>(); }
);

bool basic_dpd_engine_v7_raw_cache_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV7RawCache",
    [](){ return std::make_shared<BasicDPDEngineV7Raw<true>>(); }
);