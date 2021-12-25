#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw.hpp"

bool basic_dpd_engine_v5_f22_raw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v5_f22_raw",
    [](){ return std::make_shared<BasicDPDEngineV5F22Raw>(); }
);