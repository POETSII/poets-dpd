#if defined(__x86_64__)

#include "dpd/engines/simd/dpd_engine_avx2_half_merge_tbb.hpp"

bool dpd_engine_avx2_half_merge_tbb_hash = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_hash",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<DPDEngineAVX2HalfMergeTBBConfigBase>>(); }
);

#endif