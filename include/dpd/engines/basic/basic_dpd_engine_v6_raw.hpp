#ifndef basic_dpd_engine_v6_raw_hpp
#define basic_dpd_engine_v6_raw_hpp

#include "dpd/engines/basic/basic_dpd_engine.hpp"

#include "dpd/engines/basic/basic_dpd_engine_v6_raw_handlers.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>
#include <variant>
#include <random>

/*
    The main difference here is that the devices are intended to stay active
    between calls to Run. This is to avoid needing to reboot the tinsel
    devices, as it is too slow.

    All devices start with no resident beads. These are injected in from the
    host first before execution starts. Once all beads are loaded, a single
    message is injected, which is then fanned out within the system. It then
    moves passed the barrier, and proceeds as normal. When it finishes,
    the resident beads "leave" the cell as they are output, and it hits
    the barrier with no beads. 
*/
class BasicDPDEngineV6Raw
    : public BasicDPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=12;

    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;

    using Handlers = BasicDPDEngineV6RawHandlers;

    using OutputFlags = Handlers::OutputFlags;
    using raw_bead_view_t = Handlers::raw_bead_view_t;
    using raw_angle_bond_info_t = Handlers::raw_angle_bond_info_t;
    using raw_bead_resident_t = Handlers::raw_bead_resident_t;
    using raw_cached_bond_t = Handlers::raw_cached_bond_t;
    using raw_force_input_t = Handlers::raw_force_input_t;
    using raw_begin_t = Handlers::raw_begin_t;
    using device_state_t = Handlers::device_state_t;

    struct message_t
    {
        device_state_t *dst;
        Handlers::OutputFlags pin;
        std::variant<std::monostate,raw_bead_resident_t,raw_bead_view_t,raw_force_input_t,raw_begin_t> payload;
    };
    
    std::string CanSupport(const WorldState *s) const override
    {
        if(s->bead_types.size()>1){
            double diss=s->interactions.at(1).dissipative;
            for(unsigned i=0; i<s->bead_types.size(); i++){
                for(unsigned j=0; j<s->bead_types.size(); j++){
                    if( s->interactions[i*s->bead_types.size()+j].dissipative!=diss){
                        return "Dissipative strength must be uniform";
                    }
                }
            }
        }

        return BasicDPDEngine::CanSupport(s);
    }

    std::vector<device_state_t> m_devices;
    std::unordered_map<vec3i_t,device_state_t*> m_location_to_device;
    std::unordered_map<device_state_t*,std::vector<device_state_t*>> m_neighbour_map;
    std::unordered_map<device_state_t*,std::vector<device_state_t*>> m_begin_fanout_map;
    std::unordered_map<uint32_t,uint32_t> m_bead_hash_to_original_id;

    void set_bead_hash(raw_bead_resident_t &b, bool is_monomer, unsigned polymer_id, unsigned polymer_offset, unsigned bead_type)
    {
        b.id=bead_hash_construct(bead_type, is_monomer, polymer_id, polymer_offset);
    }

    void Attach(WorldState *state) override
    {
        m_devices.clear();
        m_location_to_device.clear();
        m_neighbour_map.clear();
        m_begin_fanout_map.clear();
        m_bead_hash_to_original_id.clear();

        BasicDPDEngine::Attach(state);

        if(state){
            m_location_to_device.reserve(m_cells.size());
            m_neighbour_map.reserve(m_cells.size());
            m_begin_fanout_map.reserve(m_cells.size());

            m_devices.resize(m_cells.size());
            for(unsigned i=0; i<m_cells.size(); i++){
                device_state_t &dst=m_devices[i];
                const cell_t &src=m_cells[i];

                m_box.extract(dst.box);
                dst.dt=m_state->dt;
                dst.inv_root_dt=recip_pow_half(m_state->dt);
                dst.bond_r0=m_bond_r0;
                dst.bond_kappa=m_bond_kappa;
                for(unsigned i=0; i<m_state->bead_types.size(); i++){
                    for(unsigned j=0; j<m_state->bead_types.size(); j++){
                        dst.conservative[i*MAX_BEAD_TYPES+j]=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                    }
                }
                dst.sqrt_dissipative=pow_half(m_state->interactions[0].dissipative);
                dst.t_hash = m_t_hash;
                dst.t_seed = m_state->seed;
                src.location.extract(dst.location);
                m_location_to_device[vec3i_t(src.location)]=&dst;

                auto &nhood = m_neighbour_map[&dst];
                nhood.reserve(src.neighbours.size());
                for(unsigned neighbour_index : src.neighbours){
                    nhood.push_back( &m_devices[neighbour_index] );
                }

                // Build a binary tree routed at 0
                unsigned parent=i/2;
                if(i!=0){
                    m_begin_fanout_map[&m_devices[parent]].push_back(&dst);
                }
            }

            for(const auto &b : state->beads){
                m_bead_hash_to_original_id[ b.get_hash_code() ] = b.bead_id;
            }
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        import_beads();

        for(auto &device : m_devices){
            device.phase=Handlers::Idle;
            device.steps_todo=0;
            device.share_todo=0;
            device.t_seed=m_state->seed;
            // Need to walk it forwards by one step
            device.t_hash=get_t_hash(device.t, device.t_seed);
            device.t=m_state->t;
            make_bag_wrapper(device.resident).clear();
            make_bag_wrapper(device.force_outgoing).clear();
            make_bag_wrapper(device.migrate_outgoing).clear();
        }

        std::vector<message_t> messages;
        for(const auto &b : m_state->beads){
            raw_bead_resident_t bb;
            import_bead(bb, b);

            // Pre-correct one-step forwards in space
            dpd_maths_core_half_step_raw::update_pos((float)m_state->dt, m_box, bb);

            message_t msg;

            vec3i_t loc=floor(bb.x);
            msg.dst = m_location_to_device[loc];
            msg.pin = Handlers::OutputFlags::RTS_INDEX_output; // For convenience we overload input as output
            msg.payload = bb;

            messages.push_back(msg);
        }

        {
            message_t msg;
            msg.dst=m_location_to_device[{0,0,0}];
            msg.pin=Handlers::OutputFlags::RTS_INDEX_begin;
            raw_begin_t payload;
            payload.num_steps = nSteps;
            msg.payload=payload;

            messages.push_back(msg);
        }

        std::vector<raw_bead_resident_t> outputs=step_all(m_devices, m_neighbour_map, m_begin_fanout_map, messages);

        for(auto &output : outputs){
            // post correct one step forwards in mom
            dpd_maths_core_half_step_raw::update_mom<float,raw_bead_resident_t>((float)m_state->dt, output);
            
            auto hash=Handlers::get_hash_code(output.id);
            auto bead_id=m_bead_hash_to_original_id.at(hash);
            auto &dst=m_state->beads.at(bead_id);
            assert(dst.get_hash_code() == hash);

            dst.x.assign(output.x);
            dst.v.assign(output.v);
            dst.f.assign(output.f);
        }

        for(unsigned i=0; i<nSteps; i++){
            m_state->t += 1;
        }
    }

    virtual std::vector<raw_bead_resident_t> step_all(
        std::vector<device_state_t> &states,
        std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map,
        std::unordered_map<device_state_t*,std::vector<device_state_t*>> &begin_fanout_map,
        std::vector<message_t> &messages
    ) {
        std::vector<message_t> messages_next;
        std::vector<raw_bead_resident_t> outputs;

        std::mt19937_64 rng;

        while(1){
            for(const auto &message : messages){
                if(rng()&1){
                    messages_next.push_back(message);
                }else{
                    switch(message.pin){
                    case Handlers::OutputFlags::RTS_INDEX_begin:
                        assert(message.dst);
                        Handlers::on_recv_begin(*message.dst, std::get<raw_begin_t>(message.payload));
                        break;
                    
                    case Handlers::OutputFlags::RTS_INDEX_output:
                        if(message.dst){
                            // This is input-from-the-host, which uses the output pin
                            Handlers::on_recv_input(*message.dst, std::get<raw_bead_resident_t>(message.payload));
                        }else{
                            // proper output to host
                            outputs.push_back(std::get<raw_bead_resident_t>(message.payload));
                            if(outputs.size() > m_state->beads.size()){
                                throw std::logic_error("More outputs than beads.");
                            }
                        }
                        break;

                    case Handlers::OutputFlags::RTS_INDEX_migrate:
                        assert(message.dst);
                        Handlers::on_recv_migrate(*message.dst, std::get<raw_bead_resident_t>(message.payload));
                        break;

                    case Handlers::OutputFlags::RTS_INDEX_force:
                        assert(message.dst);
                        Handlers::on_recv_force(*message.dst, std::get<raw_force_input_t>(message.payload));
                        break;

                    case Handlers::OutputFlags::RTS_INDEX_share:
                        assert(message.dst);
                        Handlers::on_recv_share(*message.dst, std::get<raw_bead_view_t>(message.payload));
                        break;

                    default:
                        throw std::logic_error("Unknown message type.");
                    }
                }
            }
            messages.clear();
            std::swap(messages, messages_next);
            bool active=!messages.empty();

            for(auto &state : states){
                auto rts=Handlers::calc_rts(state);
                if(rts==0){
                    continue;
                }
                active=true;
                if(rng()&1){
                    continue;
                }
                decltype(message_t::payload) payload;
                Handlers::OutputFlags pin;
                bool is_local_broadcast=true; // As opposed to going to host
                bool is_begin_fanout=false; // Following the expansion tree for begin signal
                if(rts & Handlers::OutputFlags::RTS_FLAG_force){
                    raw_force_input_t tmp;
                    Handlers::on_send_force(state, tmp);
                    pin=Handlers::OutputFlags::RTS_INDEX_force;
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_migrate){
                    raw_bead_resident_t tmp;
                    Handlers::on_send_migrate(state, tmp);
                    pin=Handlers::OutputFlags::RTS_INDEX_migrate;
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_share){
                    raw_bead_view_t tmp;
                    Handlers::on_send_share(state, tmp);
                    pin=Handlers::OutputFlags::RTS_INDEX_share;
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_begin){
                    raw_begin_t tmp;
                    Handlers::on_send_begin(state, tmp);
                    pin=Handlers::OutputFlags::RTS_INDEX_begin;
                    payload=tmp;
                    is_local_broadcast=false;
                    is_begin_fanout=true;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_output){
                    raw_bead_resident_t tmp;
                    Handlers::on_send_output(state, tmp);
                    pin=Handlers::OutputFlags::RTS_INDEX_output;
                    payload=tmp;
                    is_local_broadcast=false;
                }else{
                    throw std::logic_error("Unknown output flag.");
                }

                if(is_local_broadcast){
                    for(auto neighbour : neighbour_map[&state]){
                        messages.push_back({neighbour, pin, payload});
                    }
                }else if(is_begin_fanout){
                    for(auto neighbour : begin_fanout_map[&state]){
                        messages.push_back({neighbour, pin, payload});
                    }
                }else{
                    messages.push_back({nullptr, pin, payload});
                }

                rts=Handlers::calc_rts(state);
                if(rts){
                    active=true;
                }
            }

            if(active){
                continue;
            }

            if(outputs.size()==m_state->beads.size()){
                break;
            }

            // Implement barrier logic
            assert(!active);
            assert(messages.empty());

            unsigned active_yes=0;
            for(auto &state : states){
                active_yes += Handlers::on_barrier(state);
            }
            if(active_yes==0){
                break;
            }
            if(active_yes == states.size()){
                continue;
            }
            throw std::logic_error("Devices do no agree on activeness at barrier.");
        }

        for(auto &s : states){
            assert(s.phase==Handlers::Phase::Idle);
        }

        assert(outputs.size()==m_state->beads.size());

        return outputs;
    }


};

#endif
