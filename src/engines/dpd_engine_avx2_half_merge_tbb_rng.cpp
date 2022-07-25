#include "dpd/engines/simd/dpd_engine_avx2_half_merge_tbb.hpp"

struct config : DPDEngineAVX2HalfMergeTBBConfigBase
{
    static const typename dpd_maths_core_simd::Flags MathsCoreFlags = dpd_maths_core_simd::Flag_RngXorShift64Add;
};

bool dpd_engine_avx2_half_merge_tbb_rng = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<config>>(); }
);
