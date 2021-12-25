#include "dpd/engines/basic/basic_dpd_engine_v2.hpp"

bool basic_dpd_engine_v2_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v2",
    [](){ return std::make_shared<BasicDPDEngineV2>(); }
);