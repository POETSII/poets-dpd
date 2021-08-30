#ifndef basic_dpd_engine_v7_raw_tinsel_hpp
#define basic_dpd_engine_v7_raw_tinsel_hpp

#include "dpd/engines/basic/basic_dpd_engine_v7_raw.hpp"

#ifndef TINSEL
#include "POLiteSWSim_PGraph.h"
#include <iostream>
#include <unistd.h>
#endif

#include "POLiteHW.h"

template<class Impl = POLiteHW<>>
class BasicDPDEngineV7RawTinsel
    : public BasicDPDEngineV7Raw
{
public:
    using Handlers = BasicDPDEngineV7RawHandlers<false>;

    using None = typename Impl::None;

    static constexpr auto NHoodPin = Impl::Pin(0);

    struct Message
    {
        OutputFlags type;
        union{
            raw_bead_resident_t bead_resident;
            raw_bead_share_t bead_share;
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
            auto rts=Handlers::calc_rts(s->state);
            if(rts&OutputFlags::RTS_FLAG_output){
                *readyToSend = Impl::HostPin;
            }else if(rts){
                *readyToSend = NHoodPin;
            }else{
                *readyToSend = Impl::No;
            }
        }

        inline void init() {
            Handlers::on_init(s->state);
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
                Handlers::on_send_share(s->state, (raw_bead_share_t&)msg->bead_share);
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
                Handlers::on_recv_share(s->state, msg->bead_share);

            }else if(msg->type==OutputFlags::RTS_INDEX_force){
                Handlers::on_recv_force(s->state, msg->force_input);

            }else if(msg->type==OutputFlags::RTS_INDEX_migrate){
                Handlers::on_recv_migrate(s->state, msg->bead_resident);
            
            }else{
                assert(false);
            }

            update_rts();
        }

        inline bool finish(BasicDPDEngineV7RawTinsel::Message volatile*)
        { return false; }
    };

    using Thread = PThread<
          Device,
          State,     // State
          None,         // Edge label
          Message,    // Message
          1
        >;

#ifndef TINSEL
    std::shared_ptr<typename Impl::HostLink> m_hostlink;
    std::shared_ptr<typename Impl::template PGraph<Device, State, None, Message>> m_graph;
    int m_meshLenX;
    int m_meshLenY;

    BasicDPDEngineV7RawTinsel()
    {
        std::cerr<<"Construct\n";
        m_meshLenX=Impl::BoxMeshXLen;
        m_meshLenY=Impl::BoxMeshYLen;
    }

    void ensure_hostlink()
    {
        if(!m_hostlink){
            //std::cerr<<"Opening hostlink\n";
            typename Impl::HostLinkParams params;
            params.numBoxesX=m_meshLenX;
            params.numBoxesY=m_meshLenY;
            //std::cerr<<"  boxx="<<params.numBoxesX<<", boxy="<<params.numBoxesY<<"\n";
            m_hostlink = std::make_shared<typename Impl::HostLink>(params);
        }
    }

    void Attach(WorldState *state) override
    {
        ensure_hostlink();
        
        BasicDPDEngineV7Raw::Attach(state);

        //std::cerr<<"Building graph\n";

        m_graph=std::make_shared<typename Impl::template PGraph<Device, State, None, Message>>(m_meshLenX, m_meshLenY);
        auto &graph=*m_graph;

        unsigned volume=state->box[0]*state->box[1]*state->box[2];
        
        if(volume > 40*40*40){
            // TOD: Link to the number of boxes
            graph.mapVerticesToDRAM=true;
        }

        for(unsigned i=0; i<m_devices.size(); i++){
            auto id=graph.newDevice();
            if(id!=i){
                throw std::logic_error("Device ids not sequentials.");
            }
        }
        // PGraph ids are simply the index of the state in states
        auto get_device_id=[&](const device_state_t *s)
        {
            size_t index=s-&m_devices[0];
            assert(index<m_devices.size());
            return index;
        };

        // Build up connectivity
        // We only have two types of connectivity:
        // - HostPin : doing output (already setup)
        // - NHoodPin : sending to all neighbours including self
        for(unsigned i=0; i<m_devices.size(); i++){
            auto &s = m_devices[i];
            for(auto nb : m_neighbour_map[&s]){
                unsigned target=get_device_id(nb);
                graph.addEdge(i, /*PinId*/0, target);
            }
        }

        graph.map();

    }

    void step_all(
        std::vector<device_state_t> &states,
        std::unordered_map<device_state_t*, std::vector<device_state_t*>> &/*neighbour_map*/,
        unsigned /*interval_size*/,
        unsigned /*interval_count*/,
        std::function<bool(raw_bead_resident_t &output)> callback
    ) override
    {
        ensure_hostlink();

        unsigned nBeads=m_state->beads.size();

        auto &graph=*m_graph;

        for(unsigned i=0; i<states.size(); i++){
            graph.devices[i]->state.state = states[i];
            graph.devices[i]->state.state.t_hash = m_t_hash;
            graph.devices[i]->state.state.t_seed = m_state->seed;
        }

        //std::cerr<<"Writing graph\n";
        graph.write(m_hostlink.get());

        //std::cerr<<"Booting\n";
        m_hostlink->boot("bin/engines/basic_dpd_engine_v7_raw_tinsel_hw.riscv.code.v", "bin/engines/basic_dpd_engine_v7_raw_tinsel_hw.riscv.data.v");
        //std::cerr<<"Go\n";
        m_hostlink->go();

        auto now=[]()
        {
            timespec ts;
            clock_gettime(CLOCK_REALTIME, &ts);
            return ts.tv_sec+1e-9*ts.tv_nsec;
        };

        
        //std::cerr<<"Waiting for output\n";
        unsigned seen=0;
        double tStart=now();
        double tPrint=tStart+4;
        while(1){
            while(m_hostlink->pollStdOut(stderr));

            typename Impl::template PMessage<Message> msg;
            if(m_hostlink->canRecv()){
                ++seen;
                m_hostlink->recvMsg(&msg, sizeof(msg));
                if(!callback(msg.payload.bead_resident)){
                    break;
                }
                if(seen==1){
                    double tNow=now();
                    uint64_t nSteps = m_devices[0].interval_size * (uint64_t)m_devices[0].intervals_todo;
                    //std::cerr<<"First bead, t="<<(tNow-tStart)<<", nBeads="<<nBeads<<", nSteps="<<nSteps<<", bead*step/sec="<< (nBeads*nSteps) / (tNow-tStart)<<"\n";
                }
            }else{
                double tNow=now();
                if(tNow>tPrint){
                    std::cerr<<"  got "<<seen<<" out of "<<nBeads<<"\n";
                    tPrint=tStart+1.5*(tNow-tStart);
                }
                usleep(1);
            }
        }

        for(int i=0; i<1000; i++){
            (m_hostlink->pollStdOut(stderr));
        }

        m_hostlink.reset();
    }
#endif


};

#endif
