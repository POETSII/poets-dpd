#include "dpd/engines/basic/basic_dpd_engine_v3.hpp"

bool basic_dpd_engine_v3_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV3",
    [](){ return std::make_shared<BasicDPDEngineV3>(); }
);