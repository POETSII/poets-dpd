#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

#include "POLiteHW.h"

bool basic_dpd_engine_v5_raw_no_random_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v5_raw_no_random_tinsel_hw",
    [](){ return std::make_shared<BasicDPDEngineV5RawNoRandomTinsel<POLiteHW<>>>(); }
);
