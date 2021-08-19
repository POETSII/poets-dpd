#include "dpd/engines/basic/basic_dpd_engine.hpp"

bool basic_dpd_engine_registered = DPDEngineFactory::RegisterFactory(
    "BasicDPDEngine",
    [](){ return std::make_shared<BasicDPDEngine>(); }
);