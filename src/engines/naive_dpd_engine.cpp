#include "dpd/engines/naive/naive_dpd_engine.hpp"

bool naive_dpd_engine_registered = DPDEngineFactory::RegisterFactory(
    "naive_dpd_engine",
    [](){ return std::make_shared<NaiveDPDEngine<>>(); }
);