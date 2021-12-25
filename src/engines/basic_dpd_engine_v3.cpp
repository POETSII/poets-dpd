#include "dpd/engines/basic/basic_dpd_engine_v3.hpp"

bool basic_dpd_engine_v3_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v3",
    [](){ return std::make_shared<BasicDPDEngineV3>(); }
);