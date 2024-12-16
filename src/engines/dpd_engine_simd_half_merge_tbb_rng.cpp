#if 0 && defined(__arm64)

#include "dpd/engines/simd/dpd_engine_simd_half_merge_tbb.hpp"

struct config : DPDEngineSIMDHalfMergeTBBConfigBase
{
    static const typename dpd_maths_core_simd::Flags MathsCoreFlags = dpd_maths_core_simd::Flag_RngXorShift64Add;
};

bool dpd_engine_simd_half_merge_tbb_rng = DPDEngineFactory::RegisterFactory(
    "dpd_engine_simd_half_merge_tbb_rng",
    [](){ return std::make_shared<DPDEngineSIMDHalfMergeTBB<config>>(); }
);

struct config_auto : config
{
    using grid_partitioner = tbb::auto_partitioner;
};

bool dpd_engine_simd_half_merge_tbb_rng_partAuto = DPDEngineFactory::RegisterFactory(
    "dpd_engine_simd_half_merge_tbb_rng_partAuto",
    [](){ return std::make_shared<DPDEngineSIMDHalfMergeTBB<config_auto>>(); }
);

struct config_affinity : config
{
    using grid_partitioner = tbb::affinity_partitioner;
};

bool dpd_engine_simd_half_merge_tbb_rng_partAffinity = DPDEngineFactory::RegisterFactory(
    "dpd_engine_simd_half_merge_tbb_rng_partAffinity",
    [](){ return std::make_shared<DPDEngineSIMDHalfMergeTBB<config_affinity>>(); }
);

struct config_static : config
{
    using grid_partitioner = tbb::static_partitioner;
};

bool dpd_engine_simd_half_merge_tbb_rng_partStatic = DPDEngineFactory::RegisterFactory(
    "dpd_engine_simd_half_merge_tbb_rng_partStatic",
    [](){ return std::make_shared<DPDEngineSIMDHalfMergeTBB<config_static>>(); }
);

#endif