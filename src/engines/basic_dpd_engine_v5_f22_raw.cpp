#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw.hpp"

bool basic_dpd_engine_v5_f22_raw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV5F22Raw",
    [](){ return std::make_shared<BasicDPDEngineV5F22Raw>(); }
);