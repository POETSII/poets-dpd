#include "dpd/engines/naive/naive_dpd_engine_half_merge.hpp"

bool naive_dpd_engine_half_merge_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_half_merge",
    [](){ return std::make_shared<NaiveDPDEngineHalfMerge>(); }
);