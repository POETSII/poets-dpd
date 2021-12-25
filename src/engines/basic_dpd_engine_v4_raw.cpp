#include "dpd/engines/basic/basic_dpd_engine_v4_raw.hpp"

bool basic_dpd_engine_v4_raw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v4_raw",
    [](){ return std::make_shared<BasicDPDEngineV4Raw>(); }
);