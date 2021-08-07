#include "basic_dpd_engine_v5_raw_tinsel.hpp"

bool basic_dpd_engine_v5_raw_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV5RawTinselSWSim",
    [](){ return std::make_shared<BasicDPDEngineV5RawTinsel<POLiteSWSim<>>>(); }
);