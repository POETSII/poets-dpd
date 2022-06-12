#include "dpd/engines/simd/dpd_engine_avx2_half_merge_tbb.hpp"

bool dpd_engine_avx2_half_merge_tbb = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB>(); }
);
