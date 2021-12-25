#include "dpd/engines/basic/basic_dpd_engine_v8_raw.hpp"

bool basic_dpd_engine_v8_raw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v8_raw",
    [](){ return std::make_shared<BasicDPDEngineV8Raw>(); }
);