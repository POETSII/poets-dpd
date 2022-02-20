#ifndef gals_dpd_engine_v1_raw_hpp
#define gals_dpd_engine_v1_raw_hpp

#include "dpd/core/dpd_engine.hpp"
#include "dpd/engines/gals/gals_dpd_engine_v1_raw_handlers.hpp"
#include "dpd/core/make_nhood.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>
#include <variant>
#include <random>

#include "dpd/external/robin_hood.h"

class GALSDPDEngineV1Raw
    : public DPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = GALSDPDEngineV1Handlers::MAX_BEADS_PER_CELL;
    static constexpr size_t MAX_BEAD_TYPES=GALSDPDEngineV1Handlers::MAX_BEAD_TYPES;

    static constexpr bool EnableLogging = false;

    using Handlers = GALSDPDEngineV1Handlers;

    using OutputFlags = typename Handlers::OutputFlags;
    using bead_view_t = typename Handlers::bead_view_t;
    using bead_resident_t = typename Handlers::bead_resident_t;
    using device_state_t = typename Handlers::device_state_t;
    using message_t = typename Handlers::message_t;

    struct inflight_message_t
    {
        device_state_t *dst;
        message_t payload;
    };

    static void require(bool cond, const char *msg)
    {
        if(!cond){
            throw std::runtime_error(msg);
        }
    };
    
    std::string CanSupport(const WorldState *s) const override
    {
        if(s->bead_types.size()>1){
            if(s->bead_types.size() > MAX_BEAD_TYPES){
                return "Too many bead types.";
            }
        }

        for(auto p : s->polymers){
            auto pt = s->polymer_types.at(p.polymer_type);
            if(pt.bonds.size()>0){
                return "This engine does not support bonds (DPD forces only).";
            }
        }

        return "";
    }

    WorldState *m_state=0;
    std::vector<device_state_t> m_devices;
    robin_hood::unordered_flat_map<vec3i_t,device_state_t*> m_location_to_device;
    robin_hood::unordered_map<device_state_t*,std::vector<device_state_t*>> m_neighbour_map;
    robin_hood::unordered_flat_map<BeadHash,uint32_t> m_bead_hash_to_original_id;

    void Attach(WorldState *state) override
    {
        m_devices.clear();
        m_location_to_device.clear();
        m_neighbour_map.clear();
        m_bead_hash_to_original_id.clear();

        m_state=state;
        if(!state){
            return;
        }

        //fprintf(stderr, "Attaching\n");

        m_devices.resize(m_state->box[0] * m_state->box[1] * m_state->box[2]);

        m_location_to_device.reserve(m_devices.size());
        m_neighbour_map.reserve(m_devices.size());

        vec3i_t box{m_state->box};

        for(unsigned i=0; i<m_devices.size(); i++){
            device_state_t &dst=m_devices[i];

            memset(&dst, 0, sizeof(dst));

            vec3i_t loc{
                int(i % box[0]),
                int( ( i / box[0] ) % box[1] ),
                int( i / (box[0]*box[1] ))
            };

            box.extract(dst.box);
            dst.dt=m_state->dt;
            dst.scaled_inv_root_dt=pow_half(24 / m_state->dt);
            for(unsigned i=0; i<m_state->bead_types.size(); i++){
                for(unsigned j=0; j<m_state->bead_types.size(); j++){
                    dst.interactions[i*MAX_BEAD_TYPES+j].conservative=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                    dst.interactions[i*MAX_BEAD_TYPES+j].sqrt_dissipative=sqrt(m_state->interactions[i*m_state->bead_types.size()+j].dissipative);
                }
            }
            dst.t_seed = m_state->seed;
            loc.extract(dst.location);
            m_location_to_device[loc]=&dst;

            for(int d=0; d<3; d++){
                if(loc[d]==0){
                    dst.edge_mask = Handlers::EdgeMask( dst.edge_mask | (Handlers::XLeft<<(2*d)) );
                }else if(loc[d]==box[d]-1){
                    dst.edge_mask = Handlers::EdgeMask( dst.edge_mask | (Handlers::XRight<<(2*d)) );
                }
            }
        }

        auto rel_nhood_locs = make_relative_nhood_full(true);
        for(unsigned i=0; i<m_devices.size(); i++){
            device_state_t &dst=m_devices[i];
            vec3i_t loc(dst.location);
            auto abs_nhood_locs=make_absolute_nhood(rel_nhood_locs, box, loc);
            auto &nhood = m_neighbour_map[&dst];
            nhood.reserve(abs_nhood_locs.size());
            for(auto n : abs_nhood_locs){
                nhood.push_back( m_location_to_device[n] );
            }
        }

        m_bead_hash_to_original_id.reserve(m_state->beads.size());
        for(const auto &b : state->beads){
            m_bead_hash_to_original_id[ b.get_hash_code() ] = b.bead_id;
        }
        assert(m_bead_hash_to_original_id.size()==m_state->beads.size());

        //fprintf(stderr, "Attached\n");
    }

    struct output_slice
    {
        unsigned time;
        unsigned num_seen;
        std::vector<bead_resident_t> beads;
        robin_hood::unordered_flat_map<BeadHash,uint32_t> *bead_hash_to_id;

        output_slice(output_slice &&) = default;
        output_slice &operator=(output_slice &&) = default;

        output_slice(unsigned _time, robin_hood::unordered_flat_map<BeadHash,uint32_t> &_bead_hash_to_id)
            : time(_time)
            , num_seen(0)
            , bead_hash_to_id(&_bead_hash_to_id)
        {
            bead_resident_t tmp;
            tmp.id=0xFFFFFFFFul;
            beads.resize(bead_hash_to_id->size(), tmp); // All start with invalid id
        }

        void add(uint32_t t, const bead_resident_t &b)
        {
            require(t==time, "Wrong time for slice.");
            auto id=bead_hash_to_id->at(BeadHash{b.id});
            auto &dst=beads.at(id);
            if(dst.id==0xFFFFFFFFul){
                dst=b;
                num_seen += 1;
            }else{
                // This is a replicate, and should be exactly the same
                require( memcmp(&b, &dst, sizeof(b))==0, "Replicate bead does not match." );
            }
        }

        bool complete() const
        { return num_seen==beads.size(); }
    };

    std::vector<output_slice> m_slices;

protected:
    void PreRun(int interval_count, unsigned interval_size)
    {
        assert(interval_count*interval_size>0);

        for(auto &device : m_devices){
            device.t=m_state->t;
            device.interval_size=interval_size;
            device.intervals_todo=interval_count;
            device.interval_offset=interval_size-1;
            device.t_hash = get_t_hash(m_state->t, m_state->seed);
            device.max_t=m_state->t + interval_size*interval_count;
            device.rts=0;

            Handlers::clear_round_header(device.rounds[0]);
            Handlers::clear_round_header(device.rounds[1]);
            auto &curr=device.rounds[device.t&1];
            curr.t=device.t;
            auto &next=device.rounds[1-(device.t&1)];
            next.t=device.t+1;

            device.rts = Handlers::RTS_FLAG_nhood;

            Handlers::invariants(device);
        }

        for(auto &b : m_state->beads){
            bead_resident_t bb;
            bb.id=b.get_hash_code().hash;
            b.x.extract(bb.x);
            b.v.extract(bb.v);
            b.f.extract(bb.f);
            // Pre-correct one-step forwards in time, as handlers will do one too few
            dpd_maths_core_half_step_raw::update_pos((float)m_state->dt, m_state->box, bb);

            vec3i_t loc=floor(bb.x);
            auto dst = m_location_to_device.at(loc);
            auto &curr=dst->rounds[dst->t&1];
            assert(curr.t==dst->t);
            
            if(curr.num_resident+1 >= Handlers::MAX_BEADS_PER_CELL){
                fprintf(stderr, "%d beads in cell.\n", curr.num_resident+1);
                throw std::runtime_error("Initial config is too dense.");
            }

            bb.t=curr.t;
            vec3_copy(bb.src, dst->location);

            Handlers::add_resident(*dst, curr, bb);

            Handlers::invariants(*dst);
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

        unsigned done=0;
        bool aborted=false;

        auto process_slice=[&](output_slice &slice) -> bool
        {
            assert(!aborted);

            auto &outputs = slice.beads;
            for(unsigned i=0; i<outputs.size(); i++){
                auto &output=outputs[i];

                dpd_maths_core_half_step_raw::update_mom<float,bead_resident_t>((float)m_state->dt, output);

                auto &dst=m_state->beads.at(i);
                assert(dst.get_hash_code() == BeadHash{output.id});

                dst.x.assign(output.x);
                dst.v.assign(output.v);
                dst.f.assign(output.f);
            }

            for(unsigned i=0; i<interval_size; i++){
                m_state->t += 1;
            }
            done += interval_size;
            interval_count -= 1;
            //fprintf(stderr, "Finished slice and moved state->t to %u\n", m_state->t);


            bool carry_on = interval_callback() && (interval_count>0);
            aborted = !carry_on;
            return carry_on;
        };

        unsigned final_slice_t=m_state->t + interval_size*interval_count;
        unsigned next_slice_t=m_state->t + interval_size; // time of the next slice to be added to slices
        int finished_slice_t=-1;
        auto &slices=m_slices;
        slices.clear();

        /*
        fprintf(stderr, "tInit=%u, interval_size=%u, interval_count=%u\n", m_state->t, interval_size, interval_count);
        fprintf(stderr, "Output times = ");
        for(unsigned i=0; i<(unsigned)interval_count; i++){
            fprintf(stderr, " %u", m_state->t+(i+1)*interval_size);
        }
        fprintf(stderr, "\n");
        */


        auto process_output=[&](uint32_t t, bead_resident_t &output) -> bool
        {
            assert(!aborted);

            require(t <= final_slice_t, "Output is from beyond slice horizon.");

            require((int)t > finished_slice_t, "Bead comes from a finished slice.");

            unsigned slice_i=0;
            while(1){
                if(slice_i==slices.size()){
                    slices.push_back(output_slice(next_slice_t, m_bead_hash_to_original_id));
                    //fprintf(stderr, "Begin slice %u at time %u\n", slice_i, slices.back().time);
                    next_slice_t += interval_size;
                }
                output_slice &s = slices.at(slice_i);
                if(t < s.time){
                    //fprintf(stderr, "  Slice %u, time=%u, size=%u\n", slice_i, s.time, s.num_seen);
                }
                require(t >= s.time, "Time does not match a slice time.");
                if(t == s.time){
                    //fprintf(stderr, "  Slice %u, time=%u, size=%u\n", slice_i, s.time, s.num_seen);
                    bool prev_comp=s.complete();
                    s.add(t, output);
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

        //fprintf(stderr, "slices=%d\n", slices.size());
        m_slices.clear();

        return done;
    }

    virtual void Run(unsigned nSteps) override
    {
        Run(1, nSteps, []() -> bool { return false; });
    }

protected:
    virtual void step_all(
        std::vector<device_state_t> &states,
        robin_hood::unordered_map<device_state_t*, std::vector<device_state_t*>> &neighbour_map,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool(uint32_t t,bead_resident_t &output)> callback
    )
    {
        //fprintf(stderr, "Stepping\n");

        assert(interval_count * interval_size > 0);

        unsigned nBeads = m_state->beads.size();

        std::vector<inflight_message_t> messages;
        std::vector<inflight_message_t> messages_next;

        std::mt19937_64 rng;

        for(auto &state : states){
            Handlers::on_init(state);
        }

        unsigned next_msg_t = interval_size;
        while(1){
            for(auto &message : messages){
                if((rng()&7)==0){
                    messages_next.push_back(message);
                }else if(message.dst==0){
                    assert(message.payload.type=Handlers::Output);
                    //fprintf(stderr, "Msg for t=%d\n", message.payload.t);
                    bool carry_on=callback(message.payload.t, message.payload.bead);
                    if(!carry_on){
                        //fprintf(stderr, "Break\n");
                        m_slices.clear();
                        return; // Quit the whole loop
                    }
                }else{
                    Handlers::on_recv(*message.dst, message.payload);
                }
            }
            messages.clear();
            std::swap(messages, messages_next);
            bool active=!messages.empty();

            for(auto &state : states){
                auto rts=Handlers::calc_rts(state);
                auto &curr=state.rounds[state.t&1];
                //std::cerr<<"  "<<state.location[0]<<", "<<state.location[1]<<", "<<state.location[2]<<" : t="<<state.t<<", max_t="<<state.max_t<<", rts="<<rts<<", num_resident="<<curr.num_resident<<", sent_exp_views="<<(int)curr.have_sent_expected_views<<", send_exp_migr="<<(int)curr.have_sent_expected_migrations<<", exp_migr_recv="<<curr.expected_migrations_received<<", exp_migr_delta="<<curr.expected_migrations_delta<<"\n";
                if(rts==0){
                    continue;
                }
                assert(rts==Handlers::RTS_FLAG_output || rts==Handlers::RTS_FLAG_nhood);
                active=true;
                if((rng()&7)==0){
                    continue;
                }
                message_t payload;
                memset(&payload, 0, sizeof(payload));
                Handlers::on_send(state, payload);
                
                if(rts&Handlers::RTS_FLAG_nhood){
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

            fprintf(stderr, "Barrier\n");
            unsigned nBeads=0, nSent=0;
            std::unordered_set<BeadHash> seen;
            for(unsigned i=0; i<m_devices.size(); i++){
                auto &d=m_devices[i];
                auto &curr=d.rounds[d.t&1];
                fprintf(stderr, "d=%u, t=%u, curr.t=%u, round_complete=%d\n", i, d.t, curr.t, Handlers::is_round_complete(d, curr));
                nBeads += curr.num_resident;
                for(unsigned j=0; j<curr.num_resident; j++){
                    if(!seen.insert(BeadHash{curr.beads[j].id}).second){
                        fprintf(stderr, "Dup\n");
                    }
                }
                assert(d.outputs_waiting==0);
            }
            fprintf(stderr, "nBeads=%u of %d, outputs=%u\n", nBeads, (int)m_state->beads.size(), nSent);
            fprintf(stderr, "nslices=%u\n", (unsigned)m_slices.size());
            for(unsigned i=0; i<m_slices.size(); i++){
                fprintf(stderr, "  slices[%u], %u of %d\n", i, m_slices[i].num_seen, (int)m_slices[i].beads.size());
            }
            assert(0);

            throw std::logic_error("Devices went idle without all outputs sent.");
        }
    }


};


#endif
