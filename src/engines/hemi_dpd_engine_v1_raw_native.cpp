#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw.hpp"

bool hemi_dpd_engine_v1_raw_native_registered = DPDEngineFactory::RegisterFactory(
    "hemi_dpd_engine_v1_raw_native",
    [](){ return std::make_shared<HemiDPDEngineV1Raw>(); }
);