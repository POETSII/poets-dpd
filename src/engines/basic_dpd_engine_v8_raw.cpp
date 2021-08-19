#include "dpd/engines/basic/basic_dpd_engine_v8_raw.hpp"

bool basic_dpd_engine_v8_raw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV8Raw",
    [](){ return std::make_shared<BasicDPDEngineV8Raw>(); }
);