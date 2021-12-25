#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

#include "dpd/core/struct_to_c.hpp"

#include <random>

int main(int argc, const char *argv[])
{
    BasicDPDEngineV5RawHandlers::device_state_t state;
    BasicDPDEngineV5RawHandlers::raw_bead_resident_t migrate_message;
    BasicDPDEngineV5RawHandlers::raw_bead_view_t share_message;
    BasicDPDEngineV5RawHandlers::raw_force_input_t force_message;

    std::cout<<"#include <stdint.h>\n";
    std::cout<<struct_to_c<BasicDPDEngineV5RawHandlers::device_state_t>("device_state_t")<<";\n";
    std::cout<<"static_assert(sizeof(device_state_t)=="<<sizeof(state)<<");\n";

    std::cout<<struct_to_c<BasicDPDEngineV5RawHandlers::raw_bead_resident_t>("migrate_message_t")<<";\n";
    std::cout<<"static_assert(sizeof(migrate_message_t)=="<<sizeof(migrate_message)<<");\n";

    std::cout<<struct_to_c<BasicDPDEngineV5RawHandlers::raw_bead_view_t>("share_message_t")<<";\n";
    std::cout<<"static_assert(sizeof(share_message_t)=="<<sizeof(share_message)<<");\n";

    std::cout<<struct_to_c<BasicDPDEngineV5RawHandlers::raw_force_input_t>("force_message_t")<<";\n";
    std::cout<<"static_assert(sizeof(force_message_t)=="<<sizeof(force_message)<<");\n";

    std::memset(&state, 0, sizeof(state));
    std::cout<<"device_state_t state="<<"\n";
    std::cout<<struct_to_c_init(state)<<"\n";
    std::cout<<";\n";
}
