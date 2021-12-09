#include "dpd/engines/basic/basic_dpd_engine_v8_raw_tinsel.hpp"

bool basic_dpd_engine_v8_raw_tinsel_hw_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngineV8RawTinselHW",
    [](){ return std::make_shared<BasicDPDEngineV8RawTinsel<POLiteHW<>>>(); }
);
