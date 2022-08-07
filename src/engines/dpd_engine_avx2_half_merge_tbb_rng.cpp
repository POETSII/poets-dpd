#include "dpd/engines/simd/dpd_engine_avx2_half_merge_tbb.hpp"

struct config : DPDEngineAVX2HalfMergeTBBConfigBase
{
    static const typename dpd_maths_core_simd::Flags MathsCoreFlags = dpd_maths_core_simd::Flag_RngXorShift64Add;
};

bool dpd_engine_avx2_half_merge_tbb_rng = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<config>>(); }
);

struct config_auto : config
{
    using grid_partitioner = tbb::auto_partitioner;
};

bool dpd_engine_avx2_half_merge_tbb_rng_partAuto = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng_partAuto",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<config_auto>>(); }
);

struct config_affinity : config
{
    using grid_partitioner = tbb::affinity_partitioner;
};

bool dpd_engine_avx2_half_merge_tbb_rng_partAffinity = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng_partAffinity",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<config_affinity>>(); }
);

struct config_static : config
{
    using grid_partitioner = tbb::static_partitioner;
};

bool dpd_engine_avx2_half_merge_tbb_rng_partStatic = DPDEngineFactory::RegisterFactory(
    "dpd_engine_avx2_half_merge_tbb_rng_partStatic",
    [](){ return std::make_shared<DPDEngineAVX2HalfMergeTBB<config_static>>(); }
);
