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
     Multiple time-steps within loop.

     TODO:
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

    virtual void Run(unsigned nSteps) override
    {

        import_beads();

        // Create shadow device states
        std::vector<device_state_t> states;
        std::unordered_map<device_state_t*,std::vector<device_state_t*>> neighbour_map;
        states.resize(m_cells.size());
        for(unsigned i=0; i<m_cells.size(); i++){
            device_state_t &dst=states[i];
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
            dst.phase=Handlers::PreMigrate;
            dst.steps_todo=nSteps;

            auto resident=make_bag_wrapper(dst.resident);
            for(const auto &b : src.beads){
                raw_bead_resident_t tmp;
                Handlers::copy_bead_resident(&tmp, &b);
                // Pre-correct one-step backwards in time, as handlers will do one too many
                dpd_maths_core_half_step_raw::update_mom((float)-m_state->dt, tmp);
                resident.push_back(tmp);
            }

            for(unsigned neighbour_index : src.neighbours){
                neighbour_map[&dst].push_back( &states[neighbour_index] );
            }
        
        }

        step_all(states, neighbour_map);

        //////////////////////////////////////////////
        // Export shadow back out
        for(unsigned i=0; i<m_cells.size(); i++){
            device_state_t &src=states[i];
            cell_t &dst=m_cells[i];

            auto resident=make_bag_wrapper(src.resident);
            dst.beads.clear();
            for(auto &src : resident){
                // Correct for the last mom update that wasn't done
                dpd_maths_core_half_step_raw::update_mom<float,raw_bead_resident_t>((float)m_state->dt, src);

                bead_resident_t tmp;
                Handlers::copy_bead_resident(&tmp, &src);
                dst.beads.push_back(tmp);
            }
        }

        for(unsigned i=0; i<nSteps; i++){
            m_state->t += m_state->dt;
            next_t_hash(m_state->seed);
        }

        export_beads();
    }

    void step_all(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map)
    {
        std::vector<message_t> messages;
        std::vector<message_t> messages_next;

        std::mt19937_64 rng;

        while(1){
            for(const auto &message : messages){
                if(rng()&1){
                    messages_next.push_back(message);
                }else if(std::holds_alternative<raw_force_input_t>(message.payload)){
                    Handlers::on_recv_force(*message.dst, std::get<raw_force_input_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_view_t>(message.payload)){
                    Handlers::on_recv_share(*message.dst, std::get<raw_bead_view_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_resident_t>(message.payload)){
                    Handlers::on_recv_migrate(*message.dst, std::get<raw_bead_resident_t>(message.payload));
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
                }else{
                    throw std::logic_error("Unknown output flag.");
                }
                for(auto neighbour : neighbour_map[&state]){
                    messages.push_back({neighbour, payload});
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
    }


};

#endif
