#ifndef basic_dpd_engine_v5_raw_tinsel_hpp
#define basic_dpd_engine_v5_raw_tinsel_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw.hpp"

#ifndef PDPD_TINSEL
#include "POLiteSWSim_PGraph.h"
#include <iostream>
#include <unistd.h>
#include "dpd/core/AsyncHostLink.hpp"
#endif

#include "POLite/PerfCounterAccumulator.h"

template<bool NoBonds, class Impl>
class BasicDPDEngineV5RawTinselImpl
    : public BasicDPDEngineV5RawImpl<NoBonds>
{
public:
    using Handlers = BasicDPDEngineV5RawHandlersImpl<NoBonds>;
    using Base = BasicDPDEngineV5RawImpl<NoBonds>;

    using None = typename Impl::None;

    static constexpr auto NHoodPin = Impl::Pin(0);

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
    static constexpr bool ENABLE_CORE_PERF_COUNTERS = UseDevicePerfCounters;
    static constexpr bool ENABLE_THREAD_PERF_COUNTERS = UseDevicePerfCounters;


    struct Message
    {
        typename Base::OutputFlags type;
        union{
            typename Base::raw_bead_resident_t bead_resident;
            typename Base::raw_bead_view_t bead_view;
            typename Base::raw_force_input_t force_input;
        };
    };

    struct State
    {
        typename Base::device_state_t state;
        uint32_t device_performance_counters[UseDevicePerfCounters ? NumPerfCounters : 0];
    };

    

    struct Device : Impl::template PDevice<State, None, Message> {
        using Impl::template PDevice<State, None, Message>::s;
        using Impl::template PDevice<State, None, Message>::readyToSend;

        void update_rts()
        {
            auto rts=BasicDPDEngineV5RawHandlers::calc_rts(s->state);
            if(rts&Base::OutputFlags::RTS_FLAG_output){
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
            if constexpr(UseDevicePerfCounters){
                for(unsigned i=0; i< std::size(s->device_performance_counters); i++){
                    s->device_performance_counters[i]=0;
                }
                s->device_performance_counters[LocX] = s->state.location[0];
                s->device_performance_counters[LocY] = s->state.location[1];
                s->device_performance_counters[LocZ] = s->state.location[2];
            }
            update_rts();
        }

        inline bool step() {
            if constexpr(UseDevicePerfCounters){
                if(s->state.phase!=Handlers::Migrating){
                    s->device_performance_counters[ResidentBeadCount] += s->state.resident.n;
                    uint32_t bonds=0, angles=0;
                    for(unsigned i=0; i<s->state.resident.n; i++){
                        const auto &bead=s->state.resident.elements[i];
                        if(!BeadHash{bead.id}.is_monomer()){
                            for(unsigned j=0; j<Base::MAX_BONDS_PER_BEAD; j++){
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

            if(rts&Base::OutputFlags::RTS_FLAG_share){
                Handlers::on_send_share(s->state, (typename Base::raw_bead_view_t&)msg->bead_view);
                msg->type=Base::OutputFlags::RTS_INDEX_share;

            }else if(rts&Base::Handlers::RTS_FLAG_force){
                Handlers::on_send_force(s->state, (typename Base::raw_force_input_t&)msg->force_input);
                msg->type=Base::OutputFlags::RTS_INDEX_force;

            }else if(rts&Base::OutputFlags::RTS_FLAG_migrate){
                Handlers::on_send_migrate(s->state, (typename Base::raw_bead_resident_t&)msg->bead_resident);
                msg->type=Base::OutputFlags::RTS_INDEX_migrate;

            }else if(rts&Base::OutputFlags::RTS_FLAG_output){
                Handlers::on_send_output(s->state, (typename Base::raw_bead_resident_t&)msg->bead_resident);
                //printf("Out %x\n", ((raw_bead_resident_t&)msg->bead_resident).id);
                msg->type=Base::OutputFlags::RTS_INDEX_output;

            }else{
                assert(false);
            }

            update_rts();
        }

        inline void recv(Message* msg, None* /*edge*/) {
            if(msg->type==Base::OutputFlags::RTS_INDEX_share){
                Handlers::template on_recv_share<false>(s->state, msg->bead_view);

            }else if(msg->type==Base::OutputFlags::RTS_INDEX_force){
                Handlers::on_recv_force(s->state, msg->force_input);

            }else if(msg->type==Base::OutputFlags::RTS_INDEX_migrate){
                Handlers::on_recv_migrate(s->state, msg->bead_resident);
            
            }else{
                assert(false);
            }

            update_rts();
        }

        inline bool finish(Message volatile*)
        { return false; }
    };

    using Thread = typename Impl::template PThread<
          Device,
          State,     // State
          None,         // Edge label
          Message,    // Message
          1,
          ENABLE_CORE_PERF_COUNTERS,
          ENABLE_THREAD_PERF_COUNTERS
        >;

#ifndef PDPD_TINSEL
    AsyncHostLink<Impl> m_hostlink;
    std::shared_ptr<typename Impl::template PGraph<Device, State, None, Message>> m_graph;
    int m_meshLenX;
    int m_meshLenY;
    bool m_use_device_weights;
    DPDEngine::timings_t m_timings;

    static double now()
    {
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return ts.tv_sec+1e-9*ts.tv_nsec;
    };

    BasicDPDEngineV5RawTinselImpl(bool use_device_weights=false)
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
        m_use_device_weights=use_device_weights;
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

        //std::cerr<<"Building graph\n";

        m_graph=std::make_shared<typename Impl::template PGraph<Device, State, None, Message>>(m_meshLenX, m_meshLenY);
        auto &graph=*m_graph;

        unsigned volume=state->box[0]*state->box[1]*state->box[2];
        
        if(volume > 40*40*40){
            // TODO: Link to the number of boxes
            graph.mapVerticesToDRAM=true;
        }

        for(unsigned i=0; i<Base::m_devices.size(); i++){
            auto id=graph.newDevice();
            if(id!=i){
                throw std::logic_error("Device ids not sequentials.");
            }
        }
        // PGraph ids are simply the index of the state in states
        auto get_device_id=[&](const typename Base::device_state_t *s)
        {
            size_t index=s-&Base::m_devices[0];
            assert(index<Base::m_devices.size());
            return index;
        };

        // Build up connectivity
        // We only have two types of connectivity:
        // - HostPin : doing output (already setup)
        // - NHoodPin : sending to all neighbours including self
        for(unsigned i=0; i<Base::m_devices.size(); i++){
            auto &s = Base::m_devices[i];
            for(auto nb : Base::m_neighbour_map[&s]){
                unsigned target=get_device_id(nb);
                graph.addEdge(i, /*PinId*/0, target);
            }
        }

        std::vector<unsigned> device_weights;
        if(m_use_device_weights){
            const int BASE_DEVICE_WEIGHT = 2;
            device_weights.assign(Base::m_devices.size(), BASE_DEVICE_WEIGHT);
            for(const Polymer &p : Base::m_state->polymers){
                const PolymerType &pt=Base::m_state->polymer_types[p.polymer_type];
                for(const auto &bp : pt.bond_pairs){
                    unsigned mid_off= pt.bonds.at(bp.bond_offset_head).bead_offset_tail;
                    assert(mid_off == pt.bonds.at(bp.bond_offset_tail).bead_offset_head);
                    unsigned mid_id = p.bead_ids.at(mid_off);

                    vec3i_t loc=floor(Base::m_state->beads.at(mid_id).x);
                    const auto *dst = Base::m_location_to_device[loc];
                    
                    device_weights.at(dst - &Base::m_devices[0])++;
                }
            }

            for(unsigned i=0; i<device_weights.size(); i++){
                graph.setDeviceWeight(i, device_weights[i]);
            }
        }

        graph.map();

        if(m_use_device_weights){
            std::unordered_map<unsigned,std::pair<unsigned,unsigned>> thread_sums;

            FILE *tmp=fopen("weight_mapping_cost.csv", "wt");
            fprintf(stderr, "nDevices=%u\n", (unsigned)Base::m_devices.size());

            for(unsigned i=0; i<Base::m_devices.size(); i++){
                uint32_t thread = graph.getThreadIdFromDeviceId(i);
                
                uint32_t core = (thread >> Impl::LogThreadsPerCore) << Impl::LogThreadsPerCore;
                
                unsigned weight=device_weights[i];

                auto &cs = thread_sums[thread];
                cs.first += 1;
                cs.second += weight;
            }

            for(const auto &kcs : thread_sums){
                fprintf(tmp, "%u,%u,%u,%u\n", kcs.first, (kcs.first>>Impl::LogThreadsPerCore) << Impl::LogThreadsPerCore, kcs.second.first, kcs.second.second);
            }
            fclose(tmp);
        }

        m_timings.compile = now()-tBegin;

    }

    void step_all(
        std::vector<typename Base::device_state_t> &states,
        std::unordered_map<typename Base::device_state_t*, std::vector<typename Base::device_state_t*>> &/*neighbour_map*/,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool(typename Base::raw_bead_resident_t &output)> callback
    ) override
    {
        m_hostlink.beginAquisition();

        unsigned nBeads=Base::m_state->beads.size();

        assert(m_graph);
        auto &graph=*m_graph;

        // Really this is compile
        double tBegin=now();
        for(unsigned i=0; i<states.size(); i++){
            graph.devices[i]->state.state = states[i];
            graph.devices[i]->state.state.t = Base::m_state->t;
            graph.devices[i]->state.state.t_hash = Base::m_t_hash;
            graph.devices[i]->state.state.t_seed = Base::m_state->seed;
        }
        m_timings.compile += now() - tBegin;

        tBegin = now();
        m_hostlink.completeAquisition();
        m_timings.aquire = now()-tBegin;

        tBegin = now();
        //std::cerr<<"Writing graph\n";
        graph.write(m_hostlink.get());
        //std::cerr<<"Booting\n";
        m_hostlink->boot("bin/engines/basic_dpd_engine_v5_raw_tinsel_hw.riscv.code.v", "bin/engines/basic_dpd_engine_v5_raw_tinsel_hw.riscv.data.v");
        //std::cerr<<"Go\n";
        m_timings.configure = now()-tBegin;

        PolitePerfCounterAccumulator<decltype(*m_graph),Thread> perf_counters(*m_graph, m_hostlink.get());

        m_hostlink->setStdOutFilterProc([&](uint32_t threadId, const char *line){
            return perf_counters.process_line(threadId, line);
        });

        tBegin=now();
        m_hostlink->go();
        
        //std::cerr<<"Waiting for output\n";
        unsigned seen=0;
        double tPrint=tBegin+4;
        while(1){
            while(m_hostlink->pollStdOut(stderr));

            typename Impl::template PMessage<Message> msg;
            if(m_hostlink->canRecv()){
                if(seen==0){
                    m_timings.execute_to_first_bead = now() - tBegin;
                }
                ++seen;
                m_hostlink->recvMsg(&msg, sizeof(msg));
                if(!callback(msg.payload.bead_resident)){
                    break;
                }
                
            }else{
                double tNow=now();
                if(tNow>tPrint){
                    std::cerr<<"  got "<<seen<<" out of "<<nBeads<<"\n";
                    tPrint=tBegin+1.5*(tNow-tBegin);
                }
                usleep(1); // TODO : Sigh
            }
        }
        m_timings.execute_to_last_bead=now() - tBegin;

        tBegin=now();
        fprintf(stderr, "Picking up performance counters.\n");
        while(!perf_counters.is_complete()){
            m_hostlink->pollStdOut(stderr);
            double tNow=now();
            if(tNow>tPrint){
                perf_counters.print_progress();
                tPrint=tBegin+1.5*(tNow-tBegin);
            }
        }
        if(Device::HasDevicePerfCounters){
            perf_counters.dump_combined_device_counters(
                {"ResidentBeadCount","ResidentBondCount","ResidentAngleCount","LocX","LocY","LocZ"},
                stdout
            );
        }else if(Thread::ENABLE_THREAD_PERF_COUNTERS){
            perf_counters.dump_combined_thread_counters(stdout);
        }

        while(m_hostlink->pollStdOut(stderr));
        m_timings.perf_counters = now()-tBegin;

        // TODO: this is working round the inability to restart things on a hostlink.
        m_hostlink.reset();
    }
#endif


};

template<class TImpl>
using BasicDPDEngineV5RawTinsel = BasicDPDEngineV5RawTinselImpl<false, TImpl>;

template<class TImpl>
using BasicDPDEngineV5RawNoBondsTinsel = BasicDPDEngineV5RawTinselImpl<true, TImpl>;

#endif
