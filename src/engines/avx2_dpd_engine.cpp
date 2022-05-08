#include "dpd/engines/naive/avx2_dpd_engine.hpp"

bool avx2_dpd_engine_registered = DPDEngineFactory::RegisterFactory(
    "avx2_dpd_engine",
    [](){ return std::make_shared<AVX2DPDEngine>(); }
);