#ifndef basic_dpd_engine_v5_raw_orch_hpp
#define basic_dpd_engine_v5_raw_orch_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw.hpp"

#include <iostream>
#include <unistd.h>
#include <regex>


#include "dpd/core/struct_to_c.hpp"

template<bool NoBonds>
class BasicDPDEngineV5RawOrchImpl
    : public BasicDPDEngineV5RawImpl<NoBonds>
{
private:
    using Handlers = BasicDPDEngineV5RawHandlersImpl<NoBonds>;
    using Base = BasicDPDEngineV5RawImpl<NoBonds>;

    static std::string get_shared_code()
    {
        std::string path=Handlers::THIS_HEADER;
        // TODO : less hacky way of getting graph-schema include dir
        std::string cmd="pcpp --passthru-defines --passthru-unknown-exprs -I . -I include -I ../graph_schema/include "+path;

        fprintf(stderr, "Cmd = %s\n", cmd.c_str());
        FILE *f=popen(cmd.c_str(), "r"); 
        if(!f){
            throw std::runtime_error("Problem while pre-processing with pcpp (is it installed?)");
        }

        std::string acc;
        char buffer[1024];
        while(1){
            auto done=fread(buffer, 1, sizeof(buffer), f);
            if(done==0){
                if(feof(f)){
                    break;
                }else{
                    throw std::runtime_error("Problem while pre-processing with pcpp (is it installed?)");
                }
            }

            acc.insert(acc.end(), buffer, buffer+done);
        }

        pclose(f);

        return acc;
    }

public:

    BasicDPDEngineV5RawOrchImpl()
    {
    }

    void Attach(WorldState *state) override
    {
        Base::Attach(state);



    }

    void step_all(
        std::vector<typename Handlers::device_state_t> &states,
        std::unordered_map<typename Handlers::device_state_t*, std::vector<typename Handlers::device_state_t*> > &/*neighbour_map*/,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool(typename Handlers::raw_bead_resident_t &output)> callback
    ) override
    {
        throw std::runtime_error("Orchestrator engine doesn't support live execution.");
    }

    static std::string create_graph_type_xml()
    {
        std::string header_path=Handlers::THIS_HEADER;
        std::string template_path=std::regex_replace(
            header_path,
            std::regex("_raw_handlers[.]hpp"),
            "_raw_orch.xml"
        );

        typename Handlers::device_state_t state;
        typename Handlers::raw_bead_resident_t migrate_message;
        typename Handlers::raw_bead_view_t share_message;
        typename Handlers::raw_force_input_t force_message;

        std::vector<std::pair<std::string,std::string>> replacements{
            {"__MESSAGE_TYPE_SHARE_C__",
            struct_to_c_body<typename Handlers::raw_bead_view_t>()
            },
            {"__MESSAGE_TYPE_FORCE_C__",
            struct_to_c_body<typename Handlers::raw_force_input_t>()
            },
            {"__MESSAGE_TYPE_MIGRATE_C__",
            struct_to_c_body<typename Handlers::raw_bead_resident_t>()
            },
            {"__DEVICE_STATE_C__",
            struct_to_c_body<typename Handlers::device_state_t>()
            },
            {"__SHARED_CODE_C__", get_shared_code() },
            {"__NO_BONDS__", NoBonds ? "1" : "0" },
            {"__NO_BONDS_STR__", NoBonds ? "no_bonds_" : "" }
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

        return graph_type;
    }

    void write_instance_xml(std::ostream &dst, std::string id, unsigned numTimeSteps, const std::vector<Bead> &check_beads)
    {
        unsigned intervalCount=1;
        Base::PreRun(intervalCount, numTimeSteps);

        unsigned expectedFinalTime=Base::m_state->t + numTimeSteps;
        uint32_t bead_check_sum=0;
        for(const Bead &b : Base::m_state->beads){
            bead_check_sum += b.get_hash_code().hash;
        }

        std::string graphTypeId;
        if(NoBonds){
            graphTypeId="basic_dpd_engine_v5_no_bonds_r0";
        }else{
            graphTypeId="basic_dpd_engine_v5_r0";
        }

        dst<<"<GraphInstance id='"<<id<<"' graphTypeId='"<<graphTypeId<<"'\n";
        dst<<"   P=\"{}\" >\n";
        dst<<"  <DeviceInstances>\n";

        float max_dist=0.01;
        float max_dist_squared=max_dist*max_dist;
        dst<<"    <DevI id='rr' type='reaper' P=\"{ "<<Base::m_state->beads.size()<<", "<<expectedFinalTime<<", "<<bead_check_sum<<", ";
        dst<<max_dist_squared<<", ";
        dst<<"{"<<Base::m_state->box[0]<<","<<Base::m_state->box[1]<<","<<Base::m_state->box[2]<<"},";
        dst<<check_beads.size()<<",{";
        for(unsigned i=0; i<16; i++){
            if(i!=0){
                dst<<",";
            }
            if(i < check_beads.size()){
                dst<<"{"<<check_beads[i].get_hash_code().hash<<",{"<<check_beads[i].x[0]<<","<<check_beads[i].x[1]<<","<<check_beads[i].x[2]<<"}}";
            }else{
                dst<<"{0,{0,0,0}}";
            }
        }
        dst<<"}}\" S=\"\" />\n";

        std::vector<std::string> index_to_id;
        for(unsigned i=0; i<Base::m_devices.size(); i++){
            auto &dev=Base::m_devices[i];
            auto &loc=dev.location;
            std::string id="c_"+std::to_string(loc[0])+"_"+std::to_string(loc[1])+"_"+std::to_string(loc[2]);
            index_to_id.push_back(id);

            auto state_init=struct_to_c_init(dev);

            dst<<"    <DevI id='"<<id<<"' type='cell' P=\"{}\" S=\""<<state_init<<"\" />\n";
        }

        dst<<"  </DeviceInstances>\n";
        dst<<"  <EdgeInstances>\n";

        auto get_device_id=[&](const typename Handlers::device_state_t *s)
        {
            size_t index=s-&Base::m_devices[0];
            assert(index<Base::m_devices.size());
            return index_to_id[index];
        };

        for(unsigned i=0; i<Base::m_devices.size(); i++){
            auto *device=&Base::m_devices[i];
            for(auto nb : Base::m_neighbour_map[device]){
                std::string target=get_device_id(nb);
                dst<<"   <EdgeI path=\""<<target<<":migrate_in-"<<index_to_id[i]<<":migrate_out\" />\n";
                dst<<"   <EdgeI path=\""<<target<<":share_in-"<<index_to_id[i]<<":share_out\" />\n";
                dst<<"   <EdgeI path=\""<<target<<":force_in-"<<index_to_id[i]<<":force_out\" />\n";
            }
            dst<<"    <EdgeI path=\"rr:output_in-"<<index_to_id[i]<<":output_out\" />\n";
        }

        dst<<"  </EdgeInstances>\n";
        dst<<"</GraphInstance>\n";
    }

    void write_xml(std::ostream &dst, unsigned numSteps, const std::vector<Bead> &check_beads={})
    {
        std::string graph_type=create_graph_type_xml();

        auto it=graph_type.find("</Graphs>");
        dst<<(graph_type.substr(0, it));

        write_instance_xml(dst, "blurble", numSteps, check_beads);

        dst<<"\n";
        dst<<"</Graphs>\n";
    }
};

using BasicDPDEngineV5RawOrch = BasicDPDEngineV5RawOrchImpl<false>;
using BasicDPDEngineV5RawNoBondsOrch = BasicDPDEngineV5RawOrchImpl<true>;

#endif
