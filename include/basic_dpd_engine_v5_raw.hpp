#ifndef basic_dpd_engine_v5_raw_hpp
#define basic_dpd_engine_v5_raw_hpp

#include "basic_dpd_engine.hpp"

#include "basic_dpd_engine_v5_raw_handlers.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>
#include <variant>

/*
    - Outer wrapper is true message loop.
    - Multiple time-steps within loop.
    - Output is via message.
*/
class BasicDPDEngineV5Raw
    : public BasicDPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = 8;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=8;

    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;

    using Handlers = BasicDPDEngineV5RawHandlers;

    using OutputFlags = Handlers::OutputFlags;
    using raw_bead_view_t = Handlers::raw_bead_view_t;
    using raw_angle_bond_info_t = Handlers::raw_angle_bond_info_t;
    using raw_bead_resident_t = Handlers::raw_bead_resident_t;
    using raw_cached_bond_t = Handlers::raw_cached_bond_t;
    using raw_force_input_t = Handlers::raw_force_input_t;
    using device_state_t = Handlers::device_state_t;

    struct message_t
    {
        device_state_t *dst;
        std::variant<std::monostate,raw_bead_resident_t,raw_bead_view_t,raw_force_input_t> payload;
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
    std::unordered_map<uint32_t,uint32_t> m_bead_hash_to_original_id;

    void set_bead_id(raw_bead_resident_t &b, bool is_monomer, unsigned polymer_id, unsigned polymer_offset, unsigned bead_type)
    {
        b.id=make_bead_id(is_monomer, polymer_id, polymer_offset, bead_type);
    }

    void Attach(WorldState *state) override
    {
        m_devices.clear();
        m_location_to_device.clear();
        m_neighbour_map.clear();

        BasicDPDEngine::Attach(state);

        if(state){
            m_location_to_device.reserve(m_cells.size());
            m_neighbour_map.reserve(m_cells.size());

            m_devices.resize(m_cells.size());
            for(unsigned i=0; i<m_cells.size(); i++){
                device_state_t &dst=m_devices[i];
                const cell_t &src=m_cells[i];

                m_box.extract(dst.box);
                dst.dt=m_state->dt;
                dst.inv_root_dt=1.0f/sqrtf(m_state->dt);
                dst.bond_r0=m_bond_r0;
                dst.bond_kappa=m_bond_kappa;
                for(unsigned i=0; i<m_state->bead_types.size(); i++){
                    for(unsigned j=0; j<m_state->bead_types.size(); j++){
                        dst.conservative[i*MAX_BEAD_TYPES+j]=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                    }
                }
                dst.sqrt_dissipative=sqrt(m_state->interactions[0].dissipative);
                dst.t_hash = m_t_hash;
                dst.t_seed = m_state->seed;
                src.location.extract(dst.location);
                m_location_to_device[vec3i_t(src.location)]=&dst;

                auto &nhood = m_neighbour_map[&dst];
                nhood.reserve(src.neighbours.size());
                for(unsigned neighbour_index : src.neighbours){
                    nhood.push_back( &m_devices[neighbour_index] );
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
            device.phase=Handlers::PreMigrate;
            device.steps_todo=nSteps;
            device.share_todo=0;
            device.outputs_todo=0;
            make_bag_wrapper(device.resident).clear();
            make_bag_wrapper(device.force_outgoing).clear();
            make_bag_wrapper(device.migrate_outgoing).clear();
        }

        for(const auto &b : m_state->beads){
            raw_bead_resident_t bb;
            import_bead(bb, b);

        }

        for(auto &cell : m_cells){
            for(auto &b : cell.beads){
                raw_bead_resident_t bb;
                Handlers::copy_bead_resident(&bb, &b);

                // Pre-correct one-step backwards in time, as handlers will do one too many
                dpd_maths_core_half_step_raw::update_mom((float)-m_state->dt, bb);

                vec3i_t loc=floor(b.x);
                auto &dst = m_location_to_device[loc];
                auto resident=make_bag_wrapper(dst->resident);
                resident.push_back(bb);
            }
        }

        std::vector<raw_bead_resident_t> outputs=step_all(m_devices, m_neighbour_map);

        for(auto &output : outputs){
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
            m_state->t += m_state->dt;
            next_t_hash(m_state->seed);
        }
    }

    virtual std::vector<raw_bead_resident_t> step_all(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map)
    {
        std::vector<message_t> messages;
        std::vector<message_t> messages_next;
        std::vector<raw_bead_resident_t> outputs;

        std::mt19937_64 rng;

        while(1){
            for(const auto &message : messages){
                if(rng()&1){
                    messages_next.push_back(message);
                }else if(std::holds_alternative<raw_force_input_t>(message.payload)){
                    assert(message.dst);
                    Handlers::on_recv_force(*message.dst, std::get<raw_force_input_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_view_t>(message.payload)){
                    assert(message.dst);
                    Handlers::on_recv_share(*message.dst, std::get<raw_bead_view_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_resident_t>(message.payload)){
                    if(message.dst==0){
                        outputs.push_back(std::get<raw_bead_resident_t>(message.payload));
                        if(outputs.size() > m_state->beads.size()){
                            throw std::logic_error("More outputs than beads.");
                        }
                    }else{
                        Handlers::on_recv_migrate(*message.dst, std::get<raw_bead_resident_t>(message.payload));
                    }
                }else{
                    throw std::logic_error("Unknown message type.");
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
                bool is_local_broadcast=true; // As opposed to going to host
                if(rts & Handlers::OutputFlags::RTS_FLAG_force){
                    raw_force_input_t tmp;
                    Handlers::on_send_force(state, tmp);
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_migrate){
                    raw_bead_resident_t tmp;
                    Handlers::on_send_migrate(state, tmp);
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_share){
                    raw_bead_view_t tmp;
                    Handlers::on_send_share(state, tmp);
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_output){
                    raw_bead_resident_t tmp;
                    Handlers::on_send_output(state, tmp);
                    payload=tmp;
                    is_local_broadcast=false;
                }else{
                    throw std::logic_error("Unknown output flag.");
                }
                if(is_local_broadcast){
                    for(auto neighbour : neighbour_map[&state]){
                        messages.push_back({neighbour, payload});
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

        assert(outputs.size()==m_state->beads.size());

        return outputs;
    }


};

#endif
