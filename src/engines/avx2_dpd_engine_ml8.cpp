#include "dpd/engines/naive/avx2_dpd_engine.hpp"

bool avx2_dpd_engine_ml8_registered = DPDEngineFactory::RegisterFactory(
    "avx2_dpd_engine_ml8",
    [](){ return std::make_shared<AVX2DPDEngine<8>>(); }
);
