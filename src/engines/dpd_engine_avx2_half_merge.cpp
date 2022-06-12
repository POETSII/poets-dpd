#include "dpd/engines/simd/dpd_engine_avx2_half_merge.hpp"

bool dpd_engine_avx2_half_merge = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge",
    [](){ return std::make_shared<DPDEngineAVX2HalfMerge>(); }
);
