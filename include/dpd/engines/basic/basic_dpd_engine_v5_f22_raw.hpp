#ifndef basic_dpd_engine_v5_f22_raw_hpp
#define basic_dpd_engine_v5_f22_raw_hpp

#include "dpd/engines/basic/basic_dpd_engine.hpp"

#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw_handlers.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>
#include <variant>
#include <random>

#include "dpd/external/robin_hood.h"


class BasicDPDEngineV5F22Raw
    : public BasicDPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=8;

    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;

    static constexpr bool EnableLogging = false;

    using Handlers = BasicDPDEngineV5F22RawHandlers;

    using OutputFlags = Handlers::OutputFlags;
    using raw_bead_view_f22_t = Handlers::raw_bead_view_f22_t;
    using raw_angle_bond_info_t = Handlers::raw_angle_bond_info_t;
    using raw_bead_resident_f22_t = Handlers::raw_bead_resident_f22_t;
    using raw_cached_bond_f22_t = Handlers::raw_cached_bond_f22_t;
    using raw_force_input_t = Handlers::raw_force_input_t;
    using device_state_f22_t = Handlers::device_state_f22_t;

    struct message_t
    {
        device_state_f22_t *dst;
        std::variant<std::monostate,raw_bead_resident_f22_t,raw_bead_view_f22_t,raw_force_input_t> payload;
    };
    
    std::string CanSupport(const WorldState *s) const override
    {
        if(s->bead_types.size()>1){
            if(s->bead_types.size() > MAX_BEAD_TYPES){
                return "Too many bead types.";
            }
        }

        return BasicDPDEngine::CanSupport(s);
    }

    std::vector<device_state_f22_t> m_devices;
    std::unordered_map<vec3i_t,device_state_f22_t*> m_location_to_device;
    std::unordered_map<device_state_f22_t*,std::vector<device_state_f22_t*>> m_neighbour_map;
    robin_hood::unordered_flat_map<BeadHash,uint32_t> m_bead_hash_to_original_id;

    void set_bead_id(raw_bead_resident_f22_t &b, bool is_monomer, unsigned polymer_id, unsigned polymer_offset, unsigned bead_type)
    {
        b.id=BeadHash::construct(bead_type, is_monomer, polymer_id, polymer_offset).hash;
    }

    void Attach(WorldState *state) override
    {
        m_devices.clear();
        m_location_to_device.clear();
        m_neighbour_map.clear();
        m_bead_hash_to_original_id.clear();

        BasicDPDEngine::Attach(state);

        if(state){
            m_location_to_device.reserve(m_cells.size());
            m_neighbour_map.reserve(m_cells.size());

            m_devices.resize(m_cells.size());
            for(unsigned i=0; i<m_cells.size(); i++){
                device_state_f22_t &dst=m_devices[i];
                const cell_t &src=m_cells[i];

                m_box.extract(dst.box);
                dst.dt=m_state->dt;
                dst.inv_root_dt=pow_half(24 * dpd_maths_core_half_step::kT / m_state->dt);
                dst.bond_r0=m_bond_r0;
                dst.bond_kappa=m_bond_kappa;
                for(unsigned i=0; i<m_state->bead_types.size(); i++){
                    for(unsigned j=0; j<m_state->bead_types.size(); j++){
                        dst.interactions[i*MAX_BEAD_TYPES+j].conservative=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                        dst.interactions[i*MAX_BEAD_TYPES+j].sqrt_dissipative=sqrt(m_state->interactions[i*m_state->bead_types.size()+j].dissipative);
                    }
                }
                dst.t=m_state->t;
                dst.t_hash = m_t_hash;
                dst.t_seed = m_state->seed;
                src.location.extract(dst.location_f0);
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

    struct output_slice
    {
        unsigned time;
        unsigned num_seen;
        std::vector<raw_bead_resident_f22_t> beads;
        robin_hood::unordered_flat_map<BeadHash,uint32_t> *bead_hash_to_id;

        output_slice(output_slice &&) = default;
        output_slice &operator=(output_slice &&) = default;

        output_slice(unsigned _time, robin_hood::unordered_flat_map<BeadHash,uint32_t> &_bead_hash_to_id)
            : time(_time)
            , num_seen(0)
            , bead_hash_to_id(&_bead_hash_to_id)
        {
            raw_bead_resident_f22_t tmp;
            tmp.id=0xFFFFFFFFul;
            beads.resize(bead_hash_to_id->size(), tmp); // All start with invalid id
        }

        void add(const raw_bead_resident_f22_t &b)
        {
            require(b.t==time, "Wrong time for slice.");
            auto id=bead_hash_to_id->at(Handlers::get_hash_code(b));
            auto &dst=beads.at(id);
            if(dst.id==0xFFFFFFFFul){
                dst=b;
                num_seen += 1;
            }else{
                // This is a replicate, and should be exactly the same
                require( memcmp(&b, &dst, sizeof(b))==0, "Replicate bead does not match." );
            }
            num_seen++;            
        }

        bool complete() const
        { return num_seen==beads.size(); }
    };

    virtual unsigned Run(
        int interval_count,
        unsigned interval_size,
        std::function<bool()> interval_callback
    ) {
        //std::cerr<<"Run\n";

        assert(interval_count*interval_size>0);

        import_beads();

        for(auto &device : m_devices){
            device.phase=Handlers::PreMigrate;
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

        for(auto &cell : m_cells){
            for(auto &b : cell.beads){
                raw_bead_resident_f22_t bb;
                Handlers::copy_bead_resident_plain_to_f22(&bb, &b);
                bb.t=m_state->t;

                // Pre-correct one-step backwards in time, as handlers will do one too many
                dpd_maths_core_half_step_raw::update_mom((float)-m_state->dt, bb);

                bb.checksum=calc_checksum(bb);

                vec3i_t loc{bb.x_f22[0]>>22, bb.x_f22[1]>>22, bb.x_f22[2]>>22};
                //std::cerr<<"b="<<b.x<<", bb="<<vec3i_t{bb.x_f22}<<", Loc = "<<loc<<"\n";
                auto &dst = m_location_to_device.at(loc);
                auto resident=make_bag_wrapper(dst->resident);
                resident.push_back(bb);
            }
        }

        unsigned nBeads=m_state->beads.size();

        std::vector<output_slice> slices;

        unsigned done=0;
        bool aborted=false;

        auto process_slice=[&](output_slice &slice) -> bool
        {
            assert(!aborted);

            auto &outputs = slice.beads;
            for(unsigned i=0; i<outputs.size(); i++){
                auto &output=outputs[i];

                dpd_maths_core_half_step_raw::update_mom<float,raw_bead_resident_f22_t>((float)m_state->dt, output);

                auto &dst=m_state->beads.at(i);
                assert(dst.get_hash_code() == Handlers::get_hash_code(output));

                for(int d=0; d<3; d++){
                    dst.x[d] = output.x_f22[d] * (1.0 / (1<<22));
                }
                dst.v.assign(output.v);
                dst.f.assign(output.f);
            }

            for(unsigned i=0; i<interval_size; i++){
                m_state->t += 1;
            }
            done += interval_size;
            interval_count -= 1;

            bool carry_on = interval_callback() && (interval_count>0);
            aborted = !carry_on;
            return carry_on;
        };

        unsigned final_slice_t=m_state->t + interval_size*interval_count;
        unsigned next_slice_t=m_state->t + interval_size; // time of the next slice to be added to slices
        int finished_slice_t=-1;

        auto process_output=[&](raw_bead_resident_f22_t &output) -> bool
        {
            assert(!aborted);

            require(output.t <= final_slice_t, "Output is from beyond slice horizon.");

            if((int)output.t <= finished_slice_t){
                return true; // This was a replicated output bead from a slice already finished.
            }

            unsigned slice_i=0;
            while(1){
                if(slice_i==slices.size()){
                    slices.push_back(output_slice(next_slice_t, m_bead_hash_to_original_id));
                    //fprintf(stderr, "Begin slice %u at time %u\n", slice_i, slices.back().time);
                    next_slice_t += interval_size;
                }
                output_slice &s = slices.at(slice_i);
                if(output.t < s.time){
                    fprintf(stderr, "  Bead %u, Slice %u, time=%u, size=%u\n", (unsigned)output.t, (unsigned)slice_i, (unsigned)s.time, (unsigned)s.num_seen);
                }
                //require(output.t >= s.time, "Time does not match a slice time.");
                if(output.t == s.time){
                    //fprintf(stderr, "  Slice %u, time=%u, size=%u\n", slice_i, s.time, s.num_seen);
                    bool prev_comp=s.complete();
                    s.add(output);
                    if(prev_comp!=s.complete()){
                        //fprintf(stderr, "Finished slice for time %u\n", s.time);
                    }
                    break;
                }
                ++slice_i;
            }

            while(!slices.empty() && slices.front().complete()){
                bool carry_on=process_slice(slices.front());
                finished_slice_t=slices.front().time;
                slices.erase(slices.begin());
                if(!carry_on){ // Cost should be O(nSlices), due to move 
                    return false;
                }                
            }
            return true;
        };

        step_all(m_devices, m_neighbour_map, interval_count, interval_size, process_output);

        return done;
    }

    virtual void Run(unsigned nSteps) override
    {
        Run(1, nSteps, []() -> bool { return false; });
    }

    virtual void step_all(
        std::vector<device_state_f22_t> &states,
        std::unordered_map<device_state_f22_t*, std::vector<device_state_f22_t*>> &neighbour_map,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool(raw_bead_resident_f22_t &output)> callback
    )
    {
        assert(interval_count * interval_size > 0);

        unsigned nBeads = m_state->beads.size();

        std::vector<message_t> messages;
        std::vector<message_t> messages_next;
        std::vector<raw_bead_resident_f22_t> outputs_now; // outputs for this time-step
        std::vector<raw_bead_resident_f22_t> outputs_future; // outputs for any future time-step

        std::mt19937_64 rng;

        for(auto &state : states){
            Handlers::on_init(state);
        }

        unsigned next_msg_t = interval_size;
        while(1){
            for(auto &message : messages){
                if(rng()&1){
                    messages_next.push_back(message);
                }else if(std::holds_alternative<raw_force_input_t>(message.payload)){
                    assert(message.dst);
                    Handlers::on_recv_force(*message.dst, std::get<raw_force_input_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_view_f22_t>(message.payload)){
                    assert(message.dst);
                    Handlers::on_recv_share<EnableLogging>(*message.dst, std::get<raw_bead_view_f22_t>(message.payload));
                }else if(std::holds_alternative<raw_bead_resident_f22_t>(message.payload)){
                    if(message.dst==0){
                        bool carry_on=callback(std::get<raw_bead_resident_f22_t>(message.payload));
                        if(!carry_on){
                            return; // Quit the whole loop
                        }
                    }else{
                        Handlers::on_recv_migrate(*message.dst, std::get<raw_bead_resident_f22_t>(message.payload));
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
                    raw_bead_resident_f22_t tmp;
                    Handlers::on_send_migrate(state, tmp);
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_share){
                    raw_bead_view_f22_t tmp;
                    Handlers::on_send_share(state, tmp);
                    payload=tmp;
                }else if(rts & Handlers::OutputFlags::RTS_FLAG_output){
                    raw_bead_resident_f22_t tmp;
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
    }


};

#endif
