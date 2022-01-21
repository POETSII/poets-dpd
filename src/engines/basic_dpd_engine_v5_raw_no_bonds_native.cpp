#include "dpd/engines/basic/basic_dpd_engine_v5_raw.hpp"

bool basic_dpd_engine_v5_raw_no_bonds_native_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine_v5_raw_no_bonds_native",
    [](){ return std::make_shared<BasicDPDEngineV5RawNoBonds>(); }
);