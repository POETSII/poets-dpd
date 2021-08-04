#include "naive_dpd_engine.hpp"

bool naive_dpd_engine_registered = DPDEngineFactory::RegisterFactory(
    "NaiveDPDEngine",
    [](){ return std::make_shared<NaiveDPDEngine<>>(); }
);