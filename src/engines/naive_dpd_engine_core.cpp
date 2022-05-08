#include "dpd/engines/naive/naive_dpd_engine.hpp"

bool naive_dpd_engine_core_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine_core",
    [](){ return std::make_shared<NaiveDPDEngine<true>>(); }
);