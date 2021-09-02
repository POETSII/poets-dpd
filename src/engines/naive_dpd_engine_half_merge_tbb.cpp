#include "dpd/engines/naive/naive_dpd_engine_half_merge_tbb.hpp"

bool naive_dpd_engine_half_merge_tbb_registered = DPDEngineFactory::RegisterFactory(
    "NaiveDPDEngineHalfMergeTBB",
    [](){ return std::make_shared<NaiveDPDEngineHalfMergeTBB>(); }
);