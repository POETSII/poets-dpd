#include "dpd/engines/basic/basic_dpd_engine_v6_raw_tinsel.hpp"

bool basic_dpd_engine_v6_raw_tinsel_swsim_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV6RawTinselSWSim",
    [](){ return std::make_shared<BasicDPDEngineV6RawTinsel<POLiteSWSim<>>>(); }
);