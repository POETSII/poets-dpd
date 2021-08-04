#include "basic_dpd_engine_v5_raw.hpp"

bool basic_dpd_engine_v5_raw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV5Raw",
    [](){ return std::make_shared<BasicDPDEngineV5Raw>(); }
);