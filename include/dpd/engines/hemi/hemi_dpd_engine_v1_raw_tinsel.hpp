#ifndef hemi_dpd_engine_v1_raw_tinsel_hpp
#define hemi_dpd_engine_v1_raw_tinsel_hpp

#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw.hpp"

#ifndef PDPD_TINSEL
#include "ParallelFor.h"
#include <iostream>
#include <unistd.h>
#include "dpd/core/AsyncHostLink.hpp"
#endif

template<class ImplOrig>
class HemiDPDEngineV1RawTinsel
    : public HemiDPDEngineV1Raw
{
public:
    static constexpr int NUM_PINS = 3;

    using Impl = typename ImplOrig::template rebind_num_pins<NUM_PINS>;

    using Handlers = HemiDPDEngineV1RawHandlers;
    using Base = HemiDPDEngineV1Raw;

    using None = typename Impl::None;

    static constexpr auto ForwardsPin = Impl::Pin(Handlers::RTS_INDEX_share);
    static constexpr auto BackwardsPin = Impl::Pin(Handlers::RTS_INDEX_force);
    static constexpr auto FullPin = Impl::Pin(Handlers::RTS_INDEX_migrate);

    using Message = Handlers::message_t;

    struct State
    {
        typename Base::device_state_t state;
    };

    struct Device : Impl::template PDevice<State, None, Message> {
        using Impl::template PDevice<State, None, Message>::s;
        using Impl::template PDevice<State, None, Message>::readyToSend;

        void update_rts()
        {
            auto rts=Handlers::calc_rts(s->state);
            if(rts&Base::OutputFlags::RTS_FLAG_share){
                *readyToSend =  ForwardsPin;
            }else if(rts&Base::OutputFlags::RTS_FLAG_force){
                *readyToSend =  BackwardsPin;
            }else if(rts&Base::OutputFlags::RTS_FLAG_migrate){
                //fprintf(stderr, "   rts = migrate\n");
                *readyToSend =  FullPin;
            }else if(rts){
                *readyToSend = Impl::HostPin;
            }else{
                *readyToSend = Impl::No;
            }
        }

        inline void init() {
            Handlers::on_init(s->state);
            update_rts();
            //puts("Init\n");
        }

        inline bool step() {
            //puts("Step\n");
            Base::device_state_t &state=s->state;
            bool res=Handlers::on_barrier(state);
            update_rts();
            return res;
        }

        inline void send(volatile Message* msg)
        {
            uint32_t rts=Handlers::calc_rts(s->state);
            if(rts & Handlers::RTS_FLAG_share){
                Handlers::on_send_share(s->state, (Message&)*msg);
            }else if(rts & Handlers::RTS_FLAG_force){
                Handlers::on_send_force(s->state, (Message&)*msg);
            }else if(rts & Handlers::RTS_FLAG_migrate){
                Handlers::on_send_migrate(s->state, (Message&)*msg);
            }else{
                Handlers::on_send_output(s->state, (Message&)*msg);
            }
            update_rts();
        }

        inline void recv(volatile Message* msg, None* /*edge*/) {
            if(msg->type==Handlers::RTS_INDEX_share){
                Handlers::on_recv_share(s->state, (Message&)*msg);
            }else if(msg->type==Handlers::RTS_INDEX_force){
                Handlers::on_recv_force(s->state, (Message&)*msg);
            }else if(msg->type==Handlers::RTS_INDEX_migrate){
                Handlers::on_recv_migrate(s->state, (Message&)*msg);
            }else{
                assert(0);
            }
            update_rts();
        }

        inline bool finish(Message volatile*m)
        {
            return false;
        }
    };

    using Thread = typename Impl::template PThread<
          Device,
          State,     // State
          None,         // Edge label
          Message,    // Message
          NUM_PINS    // Num Pins
        >;

#ifndef PDPD_TINSEL
    AsyncHostLink<Impl> m_hostlink;
    std::shared_ptr<typename Impl::template PGraph<Device, State, None, Message>> m_graph;
    int m_meshLenX;
    int m_meshLenY;
    DPDEngine::timings_t m_timings;

    static double now()
    {
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return ts.tv_sec+1e-9*ts.tv_nsec;
    };

    HemiDPDEngineV1RawTinsel()
    {
        std::cerr<<"Construct\n";
        m_meshLenX=Impl::BoxMeshXLen;
        m_meshLenY=Impl::BoxMeshYLen;
        if(getenv("PDPD_BOX_MESH_X")){
            int v=atoi(getenv("PDPD_BOX_MESH_X"));
            if(v<1 || v>m_meshLenX){
                throw std::runtime_error("Invalid PDPD_BOX_MESH_X");
            }
            m_meshLenX=v;
            fprintf(stderr, "Setting boxMeshLenX to %d\n", m_meshLenX);
        }
        if(getenv("PDPD_BOX_MESH_Y")){
            int v=atoi(getenv("PDPD_BOX_MESH_Y"));
            if(v<1 || v>m_meshLenY){
                throw std::runtime_error("Invalid PDPD_BOX_MESH_Y");
            }
            m_meshLenY=v;
            fprintf(stderr, "Setting boxMeshLenY to %d\n", m_meshLenY);
        }

        typename Impl::HostLinkParams params;
        params.numBoxesX=m_meshLenX;
        params.numBoxesY=m_meshLenY;
        m_hostlink.SetParams(params);
    }

    virtual bool GetTimings(DPDEngine::timings_t &timings)
    {
        timings = m_timings;
        return true;
    }

    void Attach(WorldState *state) override
    {
        if(state){
            m_hostlink.beginAquisition();
            m_timings=DPDEngine::timings_t();
        }else{
            m_hostlink.reset();
        }
        
        double tBegin=now();
        Base::Attach(state);

        if(!state){
            return;
        }

        std::cerr<<"Building graph\n";

        m_graph=std::make_shared<typename Impl::template PGraph<Device, State, None, Message>>(m_meshLenX, m_meshLenY);
        auto &graph=*m_graph;

        unsigned volume=state->box[0]*state->box[1]*state->box[2];
        
        //if(volume > 40*40*40){
            // TODO: Link to the number of boxes
            graph.mapVerticesToDRAM=true;
        //}

        auto dev0=graph.newDevices(Base::m_devices.size());
        if(dev0!=0){
            throw std::logic_error("Device ids not sequentials.");
        }

        // PGraph ids are simply the index of the state in states
        auto get_device_id=[&](const typename Base::device_state_t *s)
        {
            size_t index=s-&Base::m_devices[0];
            assert(index<Base::m_devices.size());
            return index;
        };

        // Build up connectivity
        std::cerr<<"adding edges\n";
        parallel_for_with_grain<unsigned>(0, Base::m_devices.size(), 16,
            [&](unsigned src)
            {
                auto &s = Base::m_devices[src];
                for(auto nb : Base::m_forwards_map[&s]){
                    graph.addEdgeLockedDst(src, Handlers::RTS_INDEX_share, get_device_id(nb));
                }
                for(auto nb : Base::m_backwards_map[&s]){
                    graph.addEdgeLockedDst(src, Handlers::RTS_INDEX_force, get_device_id(nb));
                }
                for(auto nb : Base::m_full_map[&s]){
                    graph.addEdgeLockedDst(src, Handlers::RTS_INDEX_migrate, get_device_id(nb));
                }
            }
        );

        std::cerr<<"Mapping\n";
        graph.map();

        m_timings.compile = now()-tBegin;

    }

    virtual void step_all(
        unsigned interval_size,
        unsigned interval_count
    ) override
    {
        m_hostlink.beginAquisition();

        unsigned nBeads=Base::m_state->beads.size();

        assert(m_graph);
        auto &graph=*m_graph;

        for(unsigned i=0; i<m_devices.size(); i++){
            graph.devices[i]->state.state = m_devices[i];
        }

        double tBegin = now();
        m_hostlink.completeAquisition();
        m_timings.aquire = now()-tBegin;

        tBegin = now();
        //std::cerr<<"Writing graph\n";
        graph.write(m_hostlink.get());
        //std::cerr<<"Booting\n";
        m_hostlink->boot("bin/engines/hemi_dpd_engine_v1_raw_tinsel_hw.riscv.code.v", "bin/engines/hemi_dpd_engine_v1_raw_tinsel_hw.riscv.data.v");
        //std::cerr<<"Go\n";
        m_timings.configure = now()-tBegin;

        tBegin=now();
        m_hostlink->go();
        
        std::cerr<<"Waiting for output\n";
        unsigned seen_outputs=0;
        unsigned seen_finish=0;
        double tPrint=tBegin+4;
        while(1){
            while(m_hostlink->pollStdOut(stderr));

            typename Impl::template PMessage<Message> msg;
            if(m_hostlink->canRecv()){
                m_hostlink->recvMsg(&msg, sizeof(msg));
                
                assert(msg.payload.type==Handlers::RTS_INDEX_output);

                if(seen_outputs==0){
                    m_timings.execute_to_first_bead = now() - tBegin;
                }

                ++seen_outputs;
                if(!m_collector.add_outputs(msg.payload.n, msg.payload.beads)){
                    break;
                }
            }else if( m_hostlink->isSim() && m_hostlink->hasTerminated() && !m_hostlink->canRecv()){
                // If we are in the simulation, the states should be synchronous with this thread,
                // and the other thread has finished

                assert(0);
            }else{
                double tNow=now();
                if(tNow>tPrint){
                    std::cerr<<"  got "<<seen_outputs<<" out of "<<nBeads<<", got "<<seen_finish<<" out of "<<m_devices.size()<<" cells\n";
                    tPrint=tBegin+1.5*(tNow-tBegin);
                }
                if(!m_hostlink->isSim()){
                    usleep(1); // TODO : Sigh
                }
            }
        }
        m_timings.execute_to_last_bead=now() - tBegin;

        tBegin=now();
        
        usleep(10000);
        while(m_hostlink->pollStdOut(stderr));
        m_timings.perf_counters = now()-tBegin;

        // TODO: this is working round the inability to restart things on a hostlink.
        m_hostlink.reset();
    }
#endif


};


#endif
