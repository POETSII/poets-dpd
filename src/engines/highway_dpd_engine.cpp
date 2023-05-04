
#include "dpd/engines/highway/dpd_engine_highway_single_thread.hpp"

bool highway_dpd_engine_hash_registered = DPDEngineFactory::RegisterFactory(
    "highway_dpd_engine_hash",
    [](){ return std::make_shared<dpd_maths_highway::HWY_NAMESPACE::DPDEngineHighwaySingleThread<true> >(); }
);

bool highway_dpd_engine_rng_registered = DPDEngineFactory::RegisterFactory(
    "highway_dpd_engine_rng",
    [](){ return std::make_shared<dpd_maths_highway::HWY_NAMESPACE::DPDEngineHighwaySingleThread<false> >(); }
);

