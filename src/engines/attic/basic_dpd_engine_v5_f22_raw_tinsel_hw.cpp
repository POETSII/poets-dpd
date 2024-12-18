#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw_tinsel.hpp"

bool basic_dpd_engine_v5_f22_raw_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v5_f22_raw_tinsel_hw",
    [](){ return std::make_shared<BasicDPDEngineV5F22RawTinsel<POLiteHW<>>>(); }
);