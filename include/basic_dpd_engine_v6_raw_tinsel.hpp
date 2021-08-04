#ifndef basic_dpd_engine_v6_raw_tinsel_hpp
#define basic_dpd_engine_v6_raw_tinsel_hpp

#include "basic_dpd_engine_v6_raw.hpp"

#ifndef TINSEL
#include "POLiteSWSim_PGraph.h"
#include <iostream>
#include <unistd.h>
#endif

#include "POLiteHW.h"

template<class Impl = POLiteHW>
class BasicDPDEngineV6RawTinsel
    : public BasicDPDEngineV6Raw
{
public:
    using Handlers = BasicDPDEngineV6RawHandlers;

    using None = typename Impl::None;

    static constexpr auto NHoodPin = Impl::Pin(0);


    struct Message
    {
        OutputFlags type;
        union{
            raw_bead_resident_t bead_resident;
            raw_bead_view_t bead_view;
            raw_force_input_t force_input;
        };
    };

    struct State
    {
        device_state_t state;
    };

    

    struct Device : Impl::template PDevice<State, None, Message> {
        using Impl::template PDevice<State, None, Message>::s;
        using Impl::template PDevice<State, None, Message>::readyToSend;

        void update_rts()
        {
            auto rts=BasicDPDEngineV6RawHandlers::calc_rts(s->state);
            if(rts&OutputFlags::RTS_FLAG_output){
                *readyToSend = Impl::HostPin;
            }else if(rts){
                *readyToSend = NHoodPin;
            }else{
                *readyToSend = Impl::No;
            }
        }

        inline void init() {
            update_rts();
        }

        inline bool step() {
            bool res= Handlers::on_barrier(s->state);
            update_rts();
            return res;
        }

        inline void send(volatile Message* msg)
        {
            uint32_t rts=Handlers::calc_rts(s->state);
            assert(rts);

            if(rts&OutputFlags::RTS_FLAG_share){
                Handlers::on_send_share(s->state, (raw_bead_view_t&)msg->bead_view);
                msg->type=OutputFlags::RTS_INDEX_share;

            }else if(rts&Handlers::RTS_FLAG_force){
                Handlers::on_send_force(s->state, (raw_force_input_t&)msg->force_input);
                msg->type=OutputFlags::RTS_INDEX_force;

            }else if(rts&OutputFlags::RTS_FLAG_migrate){
                Handlers::on_send_migrate(s->state, (raw_bead_resident_t&)msg->bead_resident);
                msg->type=OutputFlags::RTS_INDEX_migrate;

            }else if(rts&OutputFlags::RTS_FLAG_output){
                Handlers::on_send_output(s->state, (raw_bead_resident_t&)msg->bead_resident);
                msg->type=OutputFlags::RTS_INDEX_output;

            }else{
                assert(false);
            }

            update_rts();
        }

        inline void recv(Message* msg, None* /*edge*/) {
            if(msg->type==OutputFlags::RTS_INDEX_share){
                Handlers::on_recv_share(s->state, msg->bead_view);

            }else if(msg->type==OutputFlags::RTS_INDEX_force){
                Handlers::on_recv_force(s->state, msg->force_input);

            }else if(msg->type==OutputFlags::RTS_INDEX_migrate){
                Handlers::on_recv_migrate(s->state, msg->bead_resident);
            
            }else{
                assert(false);
            }

            update_rts();
        }

        inline bool finish(BasicDPDEngineV6RawTinsel::Message volatile*)
        { return false; }
    };

    using Thread = PThread<
          Device,
          State,     // State
          None,         // Edge label
          Message    // Message
        >;

#ifndef TINSEL
    std::shared_ptr<typename Impl::HostLink> m_hostlink;

    BasicDPDEngineV6RawTinsel()
    {
        std::cerr<<"Construct\n";
    }

    std::vector<raw_bead_resident_t> step_all(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map) override
    {
        if(!m_hostlink){
            std::cerr<<"Opening hostlink\n";
            m_hostlink = std::make_shared<typename Impl::HostLink>();
        }

        std::cerr<<"Building graph\n";
        typename Impl::template PGraph<Device, State, None, Message> graph;

        for(unsigned i=0; i<states.size(); i++){
            auto id=graph.newDevice();
            if(id!=i){
                throw std::logic_error("Device ids not sequentials.");
            }
        }
        // PGraph ids are simply the index of the state in states
        auto get_device_id=[&](const device_state_t *s)
        {
            size_t index=s-&states[0];
            assert(index<states.size());
            return index;
        };

        // Build up connectivity
        // We only have two types of connectivity:
        // - HostPin : doing output (already setup)
        // - NHoodPin : sending to all neighbours including self
        for(unsigned i=0; i<states.size(); i++){
            auto &s = states[i];
            for(auto nb : neighbour_map[&s]){
                unsigned target=get_device_id(nb);
                graph.addEdge(i, /*PinId*/0, target);
            }
        }

        graph.map();

        for(unsigned i=0; i<states.size(); i++){
            graph.devices[i]->state.state = states[i];
            graph.devices[i]->state.state.t_hash = m_t_hash;
            graph.devices[i]->state.state.t_seed = m_state->seed;
        }

        unsigned nBeads=m_state->beads.size();
        std::vector<raw_bead_resident_t> outputs;
        outputs.reserve(nBeads);

        std::cerr<<"Writing graph\n";
        graph.write(m_hostlink.get());

        std::cerr<<"Booting\n";
        m_hostlink->boot("bin/engines/basic_dpd_engine_v6_raw_tinsel_hw.riscv.code.v", "bin/engines/basic_dpd_engine_v6_raw_tinsel_hw.riscv.data.v");
        std::cerr<<"Go\n";
        m_hostlink->go();

        
        std::cerr<<"Waiting for output\n";
        while(outputs.size() < nBeads){
            while(m_hostlink->pollStdOut(stderr));

            typename Impl::template PMessage<Message> msg;
            if(m_hostlink->canRecv()){
                m_hostlink->recvMsg(&msg, sizeof(msg));
                outputs.push_back(msg.payload.bead_resident);
                std::cerr<<"  got "<<outputs.size()<<" out of "<<nBeads<<"\n";
            }else{
                usleep(1);
            }
        }

        m_hostlink.reset();

        return outputs;
    }
#endif


};

#endif
