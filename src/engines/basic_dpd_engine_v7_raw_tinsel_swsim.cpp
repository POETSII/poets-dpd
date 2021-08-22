#include "dpd/engines/basic/basic_dpd_engine_v7_raw_tinsel.hpp"

bool basic_dpd_engine_v7_raw_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV7RawTinselSWSim",
    [](){ return std::make_shared<BasicDPDEngineV7RawTinsel<POLiteSWSim<>>>(); }
);