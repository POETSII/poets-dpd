#ifndef basic_dpd_engine_v5_f22_raw_tinsel_hpp
#define basic_dpd_engine_v5_f22_raw_tinsel_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw.hpp"

#ifndef PDPD_TINSEL
#include "POLiteSWSim_PGraph.h"
#include <iostream>
#include <unistd.h>
#endif

#include "POLiteHW.h"

template<class Impl = POLiteHW<>>
class BasicDPDEngineV5F22RawTinsel
    : public BasicDPDEngineV5F22Raw
{
public:
    using Handlers = BasicDPDEngineV5F22RawHandlers;

    using None = typename Impl::None;

    static constexpr auto NHoodPin = Impl::Pin(0);


    struct Message
    {
        OutputFlags type;
        union{
            raw_bead_resident_f22_t bead_resident;
            raw_bead_view_f22_t bead_view;
            raw_force_input_t force_input;
        };
    };

    enum PerfCounters
    {
        ResidentBeadCount,
        ResidentBondCount,
        ResidentAngleCount,
        // Really these should not be exposed by the device, as they are static
        LocX,
        LocY,
        LocZ,

        NumPerfCounters
    };

    static constexpr bool UseDevicePerfCounters = false;
    static constexpr bool ENABLE_CORE_PERF_COUNTERS = false; //UseDevicePerfCounters;
    static constexpr bool ENABLE_THREAD_PERF_COUNTERS = false; //UseDevicePerfCounters;

    struct State
    {
        device_state_f22_t state;

        uint32_t device_performance_counters[UseDevicePerfCounters ? NumPerfCounters : 0];
    };

    

    struct Device : Impl::template PDevice<State, None, Message> {
        using Impl::template PDevice<State, None, Message>::s;
        using Impl::template PDevice<State, None, Message>::readyToSend;

        void update_rts()
        {
            auto rts=BasicDPDEngineV5F22RawHandlers::calc_rts(s->state);
            if(rts&OutputFlags::RTS_FLAG_output){
                *readyToSend = Impl::HostPin;
            }else if(rts){
                *readyToSend = NHoodPin;
            }else{
                *readyToSend = Impl::No;
            }
        }

        inline void init() {
            //printf("N %x\n", s->state.resident.elements[i].id);
            //for(int i=0; i<s->state.resident.n; i++){
            //    printf("In %x\n", s->state.resident.elements[i].id);
            //}
            Handlers::on_init(s->state);
            update_rts();

            if constexpr(UseDevicePerfCounters){
                for(unsigned i=0; i< std::size(s->device_performance_counters); i++){
                    s->device_performance_counters[i]=0;
                }
                s->device_performance_counters[LocX] = s->state.location_f0[0];
                s->device_performance_counters[LocY] = s->state.location_f0[1];
                s->device_performance_counters[LocZ] = s->state.location_f0[2];
            }
                        
        }

        inline bool step() {
            if constexpr(UseDevicePerfCounters){
                if(s->state.phase!=Handlers::Migrating){
                    s->device_performance_counters[ResidentBeadCount] += s->state.resident.n;
                    uint32_t bonds=0, angles=0;
                    for(unsigned i=0; i<s->state.resident.n; i++){
                        const auto &bead=s->state.resident.elements[i];
                        if(!BeadHash{bead.id}.is_monomer()){
                            for(unsigned j=0; j<MAX_BONDS_PER_BEAD; j++){
                                bonds += bead.bond_partners[j] != -1;
                            }
                            angles += bead.angle_bonds[0].partner_head != -1;
                        }
                    }
                    s->device_performance_counters[ResidentBondCount] += bonds;
                    s->device_performance_counters[ResidentAngleCount] += angles;
                }
            }

            bool res= Handlers::on_barrier(s->state);
            update_rts();

            return res;
        }

        inline void send(volatile Message* msg)
        {
            uint32_t rts=Handlers::calc_rts(s->state);
            assert(rts);

            if(rts&OutputFlags::RTS_FLAG_share){
                Handlers::on_send_share(s->state, (raw_bead_view_f22_t&)msg->bead_view);
                msg->type=OutputFlags::RTS_INDEX_share;

            }else if(rts&Handlers::RTS_FLAG_force){
                Handlers::on_send_force(s->state, (raw_force_input_t&)msg->force_input);
                msg->type=OutputFlags::RTS_INDEX_force;

            }else if(rts&OutputFlags::RTS_FLAG_migrate){
                Handlers::on_send_migrate(s->state, (raw_bead_resident_f22_t&)msg->bead_resident);
                msg->type=OutputFlags::RTS_INDEX_migrate;

            }else if(rts&OutputFlags::RTS_FLAG_output){
                Handlers::on_send_output(s->state, (raw_bead_resident_f22_t&)msg->bead_resident);
                //printf("Out %x\n", ((raw_bead_resident_t&)msg->bead_resident).id);
                msg->type=OutputFlags::RTS_INDEX_output;

            }else{
                assert(false);
            }

            update_rts();
        }

        inline void recv(Message* msg, None* /*edge*/) {
            
            if(msg->type==OutputFlags::RTS_INDEX_share){
                Handlers::on_recv_share<false>(s->state, msg->bead_view);

            }else if(msg->type==OutputFlags::RTS_INDEX_force){
                Handlers::on_recv_force(s->state, msg->force_input);

            }else if(msg->type==OutputFlags::RTS_INDEX_migrate){
                Handlers::on_recv_migrate(s->state, msg->bead_resident);
            
            }else{
                assert(false);
            }

            update_rts();
        }

        inline bool finish(BasicDPDEngineV5F22RawTinsel::Message volatile*)
        { return false; }
    };

    using Thread = typename Impl::template PThread<
          Device,
          State,     // State
          None,         // Edge label
          Message,    // Message
          1,   // Num pins
          ENABLE_CORE_PERF_COUNTERS,
          ENABLE_THREAD_PERF_COUNTERS
        >;

#ifndef PDPD_TINSEL
    std::shared_ptr<typename Impl::HostLink> m_hostlink;
    std::shared_ptr<typename Impl::template PGraph<Device, State, None, Message>> m_graph;
    int m_meshLenX;
    int m_meshLenY;

    BasicDPDEngineV5F22RawTinsel()
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
        
        BasicDPDEngineV5F22Raw::Attach(state);

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
        auto get_device_id=[&](const device_state_f22_t *s)
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
        std::vector<device_state_f22_t> &states,
        std::unordered_map<device_state_f22_t*, std::vector<device_state_f22_t*>> &/*neighbour_map*/,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool(raw_bead_resident_f22_t &output)> callback
    ) override
    {
        ensure_hostlink();

        unsigned nBeads=m_state->beads.size();

        auto &graph=*m_graph;

        for(unsigned i=0; i<states.size(); i++){
            graph.devices[i]->state.state = states[i];
            graph.devices[i]->state.state.t = m_state->t;
            graph.devices[i]->state.state.t_hash = m_t_hash;
            graph.devices[i]->state.state.t_seed = m_state->seed;
            if(states[i].resident.n==0){
                fprintf(stderr, "device %u at %u,%u,%u is empty\n", i, states[i].location_f0[0], states[i].location_f0[1], states[i].location_f0[2]);
            }
        }

        //std::cerr<<"Writing graph\n";
        graph.write(m_hostlink.get());

        //std::cerr<<"Booting\n";
        m_hostlink->boot("bin/engines/basic_dpd_engine_v5_f22_raw_tinsel_hw.riscv.code.v", "bin/engines/basic_dpd_engine_v5_f22_raw_tinsel_hw.riscv.data.v");
        //std::cerr<<"Go\n";
        m_hostlink->go();

        auto now=[]()
        {
            timespec ts;
            clock_gettime(CLOCK_REALTIME, &ts);
            return ts.tv_sec+1e-9*ts.tv_nsec;
        };

        /*
        std::vector<std::vector<uint32_t>> perf_counters(states.size(), std::vector<uint32_t>(Device::NumDevicePerfCounters,0));
        unsigned perf_counters_seen=0;

        auto filter=[&](uint32_t threadId, const char *line) -> bool
        {
            // Pattern is: "DPC:%x,%x,%x,%x\n",  threadId,deviceOffset,counterOffset,counterValue

            int len=strlen(line);
            if(len < 4){
                return false;
            }
            if(strncmp("DPC:",line,4)){
                return false;
            }

            unsigned gotThreadId, deviceOffset, counterOffset,counterValue;
            if(4!=sscanf(line, "DPC:%x,%x,%x,%x", &gotThreadId, &deviceOffset, &counterOffset, &counterValue)){
                return false;
            }
            if(gotThreadId!=threadId){
                fprintf(stderr, "Corrupted DPCs.\n");
                exit(1);
            }

            unsigned deviceId=graph.getDeviceId(threadId, deviceOffset);
            perf_counters.at(deviceId).at(counterOffset)=counterValue;

            perf_counters_seen++;

            return true;
        };

        m_hostlink->setStdOutFilterProc(filter);

        */

        PolitePerfCounterAccumulator<decltype(*m_graph),Thread> perf_counters(*m_graph, m_hostlink.get());

        m_hostlink->setStdOutFilterProc([&](uint32_t threadId, const char *line){
            return perf_counters.process_line(threadId, line);
        });
        
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

        fprintf(stderr, "Picking up performance counters.\n");
        while(!perf_counters.is_complete()){
            m_hostlink->pollStdOut(stderr);
            perf_counters.print_progress();
        }

        /*fprintf(stdout, "DeviceId,ThreadId,CoreId,x,y,z,SendTime,SendCount,RecvTime,RecvCount,BarrierGap,ResidentBeads,ResidentBonds,ResidentAngles\n");
        for(unsigned i=0; i<states.size(); i++){
            auto &state=states[i];
            fprintf(stdout, "%u,%u,%u,%u,%u,%u", i, graph.getThreadIdFromDeviceId(i), graph.getThreadIdFromDeviceId(i)>>16, state.location_f0[0], state.location_f0[1], state.location_f0[2]);
            for(unsigned j=0; j<perf_counters[i].size(); j++){
                fprintf(stdout, ",%u", perf_counters[i][j]);
            }
            fprintf(stdout, "\n");
        }
        */

       if(Device::HasDevicePerfCounters){
        perf_counters.dump_combined_device_counters(
            {"ResidentBeadCount","ResidentBondCount","ResidentAngleCount","LocX","LocY","LocZ"},
            stdout
        );
       }else if(Thread::ENABLE_THREAD_PERF_COUNTERS){
           perf_counters.dump_combined_thread_counters(stdout);
       }

        m_hostlink.reset();
    }
#endif


};

#endif
