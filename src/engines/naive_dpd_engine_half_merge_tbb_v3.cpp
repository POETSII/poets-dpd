#include "dpd/engines/naive/naive_dpd_engine_half_merge_tbb_v3.hpp"

bool naive_dpd_engine_half_merge_tbb_v3_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_half_merge_tbb_v3",
    [](){ return std::make_shared<NaiveDPDEngineHalfMergeTBBV3>(); }
);