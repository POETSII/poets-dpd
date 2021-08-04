#include "naive_dpd_engine.hpp"

bool naive_dpd_engine_core_registered = DPDEngineFactory::RegisterFactory(
    "NaiveDPDEngineCore",
    [](){ return std::make_shared<NaiveDPDEngine<true>>(); }
);