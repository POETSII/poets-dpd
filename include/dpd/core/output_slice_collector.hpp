#ifndef output_slice_collector_hpp
#define output_slice_collector_hpp

#include <vector>

#include "dpd/core/dpd_state.hpp"
#include "dpd/external/robin_hood.h"

#ifndef PDPD_TINSEL
#include "ParallelFor.h"
#endif


template<class TBead>
class OutputSliceCollector
{
    static void require(bool cond, const std::string &msg)
    {
        if(!cond){
            throw std::runtime_error(msg);
        }
    }

    struct output_slice
    {
        OutputSliceCollector *parent;
        unsigned time;
        std::vector<TBead> beads;
        
        output_slice(output_slice &&) = default;
        output_slice &operator=(output_slice &&) = default;

        output_slice(OutputSliceCollector *_parent, unsigned _time)
            : parent(_parent)
            , time(_time)
        {
            beads.reserve(parent->m_num_beads);
        }

        void add(const TBead &b)
        {
            assert(b.t==time);
            beads.push_back(b);
        }

        bool complete() const
        { return beads.size() == parent->m_num_beads; }
    };

    unsigned m_num_beads;
    unsigned m_interval_size;
    unsigned m_final_slice_t; // Time for the final slice  
    unsigned m_done;
    bool m_aborted;
    WorldState *m_state;
    const robin_hood::unordered_flat_map<BeadHash,uint32_t> *m_bead_hash_to_id;
    std::vector<output_slice> m_slices;

    std::function<bool()> m_callback;

    void process_output_impl(const TBead &output)
    {
        assert(!m_aborted);
        assert(!m_slices.empty());

        // Vast majority just take the first branch
        if( output.t==m_slices[0].time ){
            m_slices[0].add(output);
            return;
        }

        assert(output.t <= m_final_slice_t);


        // In "normal" cases this should capture almost all mild run-ahead
        for(unsigned i=1; i<m_slices.size(); i++){
            if(m_slices[i].time == output.t){
                m_slices[i].add(output);
                return;
            }
        }

        // Either the t is from an incorrect slice, or a new t is needed.
        // This is all possible in non-error cases, but a rare path
        require( output.t > m_slices.back().time, "Bead time is not aligned to an output slice." );
        require( output.t <= m_final_slice_t, "Bead is from after final slice.");

        while(m_slices.back().time < output.t){
            m_slices.push_back({ this, m_slices.back().time+m_interval_size });
        }

        require(m_slices.back().time == output.t, "Bead time is not aligned to an output slice.");
        
        m_slices.back().add(output);
    }

    bool process_slice(output_slice &slice);

public:
    OutputSliceCollector()
        : m_aborted(true)
        , m_state(0)
        
    {}

    OutputSliceCollector(
        WorldState *state,
        const robin_hood::unordered_flat_map<BeadHash,uint32_t> *bead_hash_to_id,
        unsigned interval_size,
        unsigned interval_count,
        std::function<bool()> callback
    )
        : m_num_beads(state->beads.size())
        , m_interval_size(interval_size)
        , m_final_slice_t( state->t + interval_size * interval_count )
        , m_done(0)
        , m_aborted(false)
        , m_state(state)
        , m_bead_hash_to_id(bead_hash_to_id)
        , m_callback(callback)
    {
        require( (interval_size * interval_count)>0, "Empty interval");

        m_slices.push_back({this, unsigned(m_state->t + m_interval_size)});
    }

    unsigned get_done() const
    {
        return m_done;
    }

    bool add_outputs(unsigned n, const TBead *beads)
    {
        assert(!m_slices.empty());

        for(unsigned i=0; i<n; i++){
            process_output_impl(beads[i]);
        }

        while(m_slices.front().complete()){
            bool carry_on=process_slice(m_slices.front());
            if(m_slices.front().time == m_final_slice_t){
                assert(m_slices.size()==1);
                carry_on = false;
            }
            if(!carry_on){
                return false;
            }

            assert(m_slices.front().time < m_final_slice_t);
            if(m_slices.size()==1){
                m_slices.push_back({this, m_slices.back().time+m_interval_size});
                assert(m_slices.back().time <= m_final_slice_t);
            }
            m_slices.erase(m_slices.begin()); // Cost should be O(nSlices), due to move, and nSlices should be small
        }

        assert(!m_slices.empty());
        return true;
    }
};

#ifndef PDPD_TINSEL
template<class TBead>
bool OutputSliceCollector<TBead>::process_slice(output_slice &slice)
{
    assert(!m_aborted);

    auto &outputs = slice.beads;
    parallel_for_with_grain<unsigned>(0, outputs.size(), 256, [&](unsigned i){
        auto &output=outputs[i];
        dpd_maths_core_half_step_raw::update_mom<float,TBead>((float)m_state->dt, output);

        unsigned index=m_bead_hash_to_id->at(BeadHash{output.id});

        auto &dst=m_state->beads.at(index);

        dst.x.assign(output.x);
        dst.v.assign(output.v);
        dst.f.assign(output.f);
    });

    m_state->t += m_interval_size;
    m_done += m_interval_size;

    bool carry_on = m_callback();
    m_aborted = !carry_on;
    return carry_on;
};
#endif

#endif