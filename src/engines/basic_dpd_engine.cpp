#include "dpd/engines/basic/basic_dpd_engine.hpp"

bool basic_dpd_engine_registered = DPDEngineFactory::RegisterFactory(
    "basic_dpd_engine",
    [](){ return std::make_shared<BasicDPDEngine>(); }
);