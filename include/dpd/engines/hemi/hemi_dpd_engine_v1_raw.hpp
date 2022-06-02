#ifndef hemi_dpd_engine_v1_raw_hpp
#define hemi_dpd_engine_v1_raw_hpp

#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw_handlers.hpp"

#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/make_nhood.hpp"

#include "dpd/core/output_slice_collector.hpp"

#include "dpd/external/robin_hood.h"

#include <iterator>
#include <unordered_map>
#include <cstring>
#include <variant>
#include <random>

class HemiDPDEngineV1Raw
    : public DPDEngine
{
public:
    using Handlers = HemiDPDEngineV1RawHandlers;

    using OutputFlags = typename Handlers::OutputFlags;
    using bead_resident_t = typename Handlers::bead_resident_t;
    using device_state_t = typename Handlers::device_state_t;
    using message_t = typename Handlers::message_t;

    static constexpr unsigned MAX_BEAD_TYPES = Handlers::MAX_BEAD_TYPES;

    struct routed_message_t
    {
        device_state_t *dst;
        message_t payload;
    };
    
    std::string CanSupport(const WorldState *s) const override
    {
        auto r=DPDEngine::CanSupport(s);
        if(!r.empty()){
            return r;
        }

        if(s->bead_types.size()>1){
            if(s->bead_types.size() > MAX_BEAD_TYPES){
                return "Too many bead types.";
            }
        }

        return {};
    }

    WorldState *m_state = 0;
    std::vector<device_state_t> m_devices;
    robin_hood::unordered_flat_map<vec3i_t,device_state_t*> m_location_to_device;
    robin_hood::unordered_map<device_state_t*,std::vector<device_state_t*>> m_forwards_map;
    robin_hood::unordered_map<device_state_t*,std::vector<device_state_t*>> m_backwards_map;
    robin_hood::unordered_map<device_state_t*,std::vector<device_state_t*>> m_full_map;
    robin_hood::unordered_flat_map<BeadHash,uint32_t> m_bead_hash_to_original_id;

    OutputSliceCollector<bead_resident_t> m_collector;

    HemiDPDEngineV1Raw()
    {}


    void Attach(WorldState *state) override
    {
        m_devices.clear();
        m_location_to_device.clear();
        m_forwards_map.clear();
        m_backwards_map.clear();
        m_full_map.clear();
        m_bead_hash_to_original_id.clear();

        m_state=state;

        if(!state){
            return;
        }

        vec3i_t box{m_state->box};

        unsigned volume=box[0]*box[1]*box[2];

        m_location_to_device.reserve(volume);
        m_forwards_map.reserve(volume);
        m_backwards_map.reserve(volume);
        m_full_map.reserve(volume);

        m_devices.resize(volume);
        unsigned index=0;
        for_each_point_in_box({0,0,0}, box, [&](vec3i_t loc){
            device_state_t &dst=m_devices[index];
            ++index;

            memzero32(dst);

            state->box.extract(dst.box);
            dst.dt=m_state->dt;
            dst.scaled_inv_root_dt=pow_half(24 / m_state->dt);
            dst.edge_bits = create_wrap_bits(&box.x[0], &loc.x[0]);
            for(unsigned i=0; i<m_state->bead_types.size(); i++){
                for(unsigned j=0; j<m_state->bead_types.size(); j++){
                    dst.interactions[i*MAX_BEAD_TYPES+j].conservative=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                    dst.interactions[i*MAX_BEAD_TYPES+j].sqrt_dissipative=sqrt(m_state->interactions[i*m_state->bead_types.size()+j].dissipative);
                }
            }
            dst.t=m_state->t;
            dst.t_hash = get_t_hash(m_state->t, m_state->seed);
            dst.t_seed = m_state->seed;
            loc.extract(dst.location);
            m_location_to_device[loc]=&dst;
        });

        auto rel_forwards=make_relative_nhood_forwards(true);
        auto rel_backwards=make_relative_nhood_inverse(rel_forwards);
        auto rel_full=make_relative_nhood_full(true);


        for_each_point_in_box({0,0,0}, box, [&](vec3i_t loc){
            device_state_t *src=m_location_to_device[loc];

            for(vec3i_t nb : make_absolute_nhood(rel_forwards, box, loc)){
                m_forwards_map[src].push_back( m_location_to_device.at(nb) );
            }
            for(vec3i_t nb : make_absolute_nhood(rel_backwards, box, loc)){
                m_backwards_map[src].push_back( m_location_to_device.at(nb) );
            }
            assert(m_full_map[src].size()==0);
            for(vec3i_t nb : make_absolute_nhood(rel_full, box, loc)){
                m_full_map[src].push_back( m_location_to_device.at(nb) );
            }
        });

        for(const auto &b : state->beads){
            m_bead_hash_to_original_id[ b.get_hash_code() ] = b.bead_id;
        }
    }

protected:
    void PreRun(int interval_count, unsigned interval_size)
    {
        assert(interval_count*interval_size>0);

        for(auto &device : m_devices){
            device.rts=0;
            device.phase=Handlers::PreMigrate;
            device.t=m_state->t;
            device.interval_size=interval_size;
            device.intervals_todo=interval_count;
            device.interval_offset=interval_size;
            device.output_reps=1;
            device.share_todo=0;
            device.outputs_todo=0;
            make_bag_wrapper(device.resident).clear();
            make_bag_wrapper(device.force_outgoing).clear();
            make_bag_wrapper(device.migrate_outgoing).clear();
        }

        for(Bead &b : m_state->beads){
            bead_resident_t bead;
            
            bead.id=b.get_hash_code().hash;
            b.x.extract(bead.x);
            for(int d=0; d<3; d++){
                assert(0.0f <= bead.x[d]);
                if(bead.x[d] == (float)m_state->box[d]){ // Handle pesky rounding due to double->float conversion
                    bead.x[d]=0.0f;
                }
                assert(bead.x[d] < m_state->box[d]);
            }

            b.v.extract(bead.v);
            b.f.extract(bead.f);
            bead.t=m_state->t;

            // Pre-correct one-step backwards in time, as handlers will do one too many
            dpd_maths_core_half_step_raw::update_mom<float,bead_resident_t>((float)-m_state->dt, bead);

            vec3i_t loc;
            vec3_floor_nn(&loc.x[0], bead.x);
            auto *dst=m_location_to_device.at(loc);
            make_bag_wrapper(dst->resident).push_back(bead);
        }
    }

public:
    virtual unsigned Run(
        int interval_count,
        unsigned interval_size,
        std::function<bool()> interval_callback
    ) {
        PreRun(interval_count, interval_size);

        unsigned nBeads=m_state->beads.size();

        m_collector=OutputSliceCollector<bead_resident_t>(
            m_state,
            &m_bead_hash_to_original_id,
            interval_size,
            interval_count,
            interval_callback
        );

        step_all(interval_count, interval_size);

        return m_collector.get_done();
    }

    virtual void Run(unsigned nSteps) override
    {
        Run(1, nSteps, []() -> bool { return false; });
    }

protected:
    virtual void step_all(
        unsigned interval_size,
        unsigned interval_count
    )
    {
        assert(interval_count * interval_size > 0);

        unsigned nBeads = m_state->beads.size();

        std::vector<routed_message_t> messages;
        std::vector<routed_message_t> messages_next;
        std::vector<bead_resident_t> outputs_now; // outputs for this time-step
        std::vector<bead_resident_t> outputs_future; // outputs for any future time-step

        std::mt19937_64 rng;

        for(auto &dev : m_devices){
            Handlers::on_init(dev);
        }

        while(1){
            for(auto &message : messages){
                if((rng()&7)==0){
                    messages_next.push_back(message);
                }else if(message.payload.type==Handlers::RTS_INDEX_force){
                    assert(message.dst);
                    Handlers::on_recv_force(*message.dst, message.payload);
                    Handlers::calc_rts(*message.dst);
                }else if(message.payload.type==Handlers::RTS_INDEX_share){
                    assert(message.dst);
                    Handlers::on_recv_share(*message.dst, message.payload);
                    Handlers::calc_rts(*message.dst);
                }else if(message.payload.type==Handlers::RTS_INDEX_output){
                    assert(message.dst==0);
                    bool carry_on=m_collector.add_outputs(message.payload.n, message.payload.beads);
                    if(!carry_on){
                        return; // Quit the whole loop
                    }
                }else if(message.payload.type==Handlers::RTS_INDEX_migrate){
                    Handlers::on_recv_migrate(*message.dst, message.payload);
                    Handlers::calc_rts(*message.dst);
                }else{
                    throw std::logic_error("Unknown message type.");
                }
            }
            messages.clear();
            std::swap(messages, messages_next);
            bool active=!messages.empty();

            for(auto &dev : m_devices){
                auto rts=Handlers::calc_rts(dev);
                if(rts==0){
                    continue;
                }
                active=true;
                if((rng()&7)==0){
                    continue;
                }
                message_t payload;
                std::vector<device_state_t*> *dests=0;
                if(rts & Handlers::OutputFlags::RTS_FLAG_force){
                    Handlers::on_send_force(dev, payload);
                    payload.type=Handlers::RTS_INDEX_force;
                    dests=&m_backwards_map[&dev];
                    Handlers::calc_rts(dev);
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_migrate){
                    Handlers::on_send_migrate(dev, payload);
                    payload.type=Handlers::RTS_INDEX_migrate;
                    dests=&m_full_map[&dev];
                    Handlers::calc_rts(dev);
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_share){
                    Handlers::on_send_share(dev, payload);
                    payload.type=Handlers::RTS_INDEX_share;
                    dests=&m_forwards_map[&dev];
                    Handlers::calc_rts(dev);
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_output){
                    Handlers::on_send_output(dev, payload);
                    payload.type=Handlers::RTS_INDEX_output;
                    Handlers::calc_rts(dev);
                }else{
                    throw std::logic_error("Unknown output flag.");
                }
                if(dests){
                    //fprintf(stderr, " ...\n");
                    for(auto neighbour : *dests){
                        messages.push_back({neighbour, payload});
                        /*if(payload.type==Handlers::RTS_INDEX_migrate){
                            fprintf(stderr, "  migrate -> [%d,%d,%d]\n", neighbour->location[0],neighbour->location[1],neighbour->location[2]);
                        }*/
                    }
                }else{
                    messages.push_back({nullptr, payload});
                }
            }

            if(active){
                continue;
            }

            // Implement barrier logic
            assert(!active);
            assert(messages.empty());

            unsigned nbeads_now=0;
            unsigned active_yes=0;
            for(auto &dev : m_devices){
                active_yes += Handlers::on_barrier(dev);
                nbeads_now += dev.resident.n + dev.migrate_outgoing.n;
                assert(nbeads_now <= nBeads);
            }
            assert(nbeads_now == nBeads);
            if(active_yes==0){
                break;
            }
            if(active_yes == m_devices.size()){
                continue;
            }
            throw std::logic_error("Devices do no agree on activeness at barrier.");
        }
    }


};


#endif
