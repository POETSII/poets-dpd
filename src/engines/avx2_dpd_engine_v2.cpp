
#if defined(__x86_64__)

#include "dpd/engines/naive/avx2_dpd_engine_v2.hpp"

bool avx2_dpd_engine_v2_registered = DPDEngineFactory::RegisterFactory(
    "avx2_dpd_engine_v2",
    [](){ return std::make_shared<AVX2DPDEngineV2>(); }
);


#endif
