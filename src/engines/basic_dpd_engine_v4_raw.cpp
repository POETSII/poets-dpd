#include "basic_dpd_engine_v4_raw.hpp"

bool basic_dpd_engine_v4_raw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV4Raw",
    [](){ return std::make_shared<BasicDPDEngineV4Raw>(); }
);