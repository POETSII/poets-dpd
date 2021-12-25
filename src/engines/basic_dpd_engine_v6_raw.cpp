#include "dpd/engines/basic/basic_dpd_engine_v6_raw.hpp"

bool basic_dpd_engine_v6_raw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v6_raw",
    [](){ return std::make_shared<BasicDPDEngineV6Raw>(); }
);