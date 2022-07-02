#include "dpd/engines/simd/dpd_engine_avx2_half_merge_tbb.hpp"

bool dpd_engine_avx2_half_merge_tbb_rng = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<
        dpd_maths_core_simd::Flags(
            dpd_maths_core_simd::Flag_RngXorShift64Add //|dpd_maths_core_simd::Flag_EnableLogging
        )
    >>(); }
);
