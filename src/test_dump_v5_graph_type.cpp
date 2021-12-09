#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

#include "dpd/core/struct_to_c.hpp"

#include <random>
#include <regex>

int main(int argc, const char *argv[])
{
    std::string header_path=BasicDPDEngineV5RawHandlers::THIS_HEADER;
    std::string template_path=std::regex_replace(
        header_path,
        std::regex("_raw_handlers[.]hpp"),
        "_raw_orch.xml"
    );

    BasicDPDEngineV5RawHandlers::device_state_t state;
    BasicDPDEngineV5RawHandlers::raw_bead_resident_t migrate_message;
    BasicDPDEngineV5RawHandlers::raw_bead_view_t share_message;
    BasicDPDEngineV5RawHandlers::raw_force_input_t force_message;

    std::vector<std::pair<std::string,std::string>> replacements{
        {"__MESSAGE_TYPE_SHARE_C__",
        struct_to_c_body<BasicDPDEngineV5RawHandlers::raw_bead_view_t>()
        },
        {"__MESSAGE_TYPE_FORCE_C__",
        struct_to_c_body<BasicDPDEngineV5RawHandlers::raw_force_input_t>()
        },
        {"__MESSAGE_TYPE_MIGRATE_C__",
        struct_to_c_body<BasicDPDEngineV5RawHandlers::raw_bead_resident_t>()
        },
        {"__DEVICE_STATE_C__",
        struct_to_c_body<BasicDPDEngineV5RawHandlers::device_state_t>()
        },
        {"__SHARED_CODE_C__", "#include \""+header_path+"\"" }
    };

    std::ifstream ifs(template_path);
    std::stringstream ss;
    ss<<ifs.rdbuf();
    std::string graph_type=ss.str();

    for(auto [name,value] : replacements){
        graph_type=std::regex_replace(
            graph_type, std::regex(name), value
        );
    }

    std::cout << graph_type;
    
}
