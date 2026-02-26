#ifndef dpd_engine_highway_single_thread_hpp
#define dpd_engine_highway_single_thread_hpp

#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/morton_codec.hpp"
#include "dpd/core/make_nhood.hpp"

#include "dpd/maths/dpd_maths_highway.hpp"
#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include <dpd/external/robin_hood.h>

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();


namespace dpd_maths_highway
{
    namespace HWY_NAMESPACE
    {
        namespace hn = hwy::HWY_NAMESPACE;

template<bool TUseHash=false>
class DPDEngineHighwaySingleThread
    : public DPDEngine
{
protected:
    static const int K = hn::MaxLanes(hn::ScalableTag<float>());

    static_assert(ExecConfig::MAX_BEAD_TYPES==8);

    struct Config_MB8 : ExecConfig
    { static constexpr bool USE_HASH = TUseHash; };  

    struct Config_MB4 : Config_MB8
    { static constexpr int MAX_BEAD_TYPES = 4; };

    struct Config_MB4_Shared : Config_MB4
    { static constexpr bool SHARED_SQRT_DISS = true; };

    struct Config_MB8_Shared : Config_MB8
    { static constexpr bool SHARED_SQRT_DISS = true; };

    struct Cell
    {
        vec3i_t pos;
        bool is_edge;
        std::vector<packed_bead> incoming; // Beads that entered cell in previous time-step
        std::array<soa_packed_neighbour,NHOOD_SIZE> nhood; // Cell 0 is self
        soa_packed_beads beads;
    };

    WorldState *m_state;

    robin_hood::unordered_flat_map<vec3i_t,uint32_t> pos_to_index;
    robin_hood::unordered_flat_map<uint32_t,uint32_t> bead_hash_to_bead_index;
    std::vector<Cell> m_cells;

    float m_rng_scale_s32_to_u;
    uint64_t m_t_hash;

    unsigned m_maxBeadTypes; // Chosen to match padding expected in maths code
    std::shared_ptr<float> m_conStrength;
    std::shared_ptr<float> m_sqrtDissStrength;
    bool m_sharedSqrtDissStrength;

    std::unique_ptr<vector_packed_beads> m_working;

    morton_codec m_morton;

    void import_beads()
    {
        for(Cell &c : m_cells){
            c.beads.clear();
        }

        float dt_2=m_state->dt*0.5f;
        for(auto b : m_state->beads){
            vec3f_t x{b.x};
            vec3i_t ix=vec3_floor(x);
            Cell &dst=m_cells[pos_to_index.at(ix)];
            packed_bead pb;
            pb.hash=b.get_hash_code().hash;
            x.extract(pb.x);
            b.v.extract(pb.v);
            b.f.extract(pb.f);

            // Correct mom backwards
            for(int d=0; d<3; d++){
                pb.v[d] -= dt_2 * pb.f[d];
            }
            dst.beads.insert(pb);

            bead_hash_to_bead_index[pb.hash]=b.bead_id;
        }
    }

    void flush_incoming(Cell &cell)
    {
        for(const auto &b : cell.incoming){
            for(int d=0; d<3; d++){
                assert(cell.pos[d] <= b.x[d]);
                assert(b.x[d] < 1+cell.pos[d]);
            }
            // Usually only one incoming, so don't need to optimise for more than one
            cell.beads.insert(b);
        }
        cell.incoming.clear();
    }

    void export_beads()
    {
        double dt_2 = m_state->dt *0.5f;

        for(Cell &cell : m_cells){
            // Could iterate over incoming without packing, but might as well do it now
            flush_incoming(cell);

            auto &beads = cell.beads;
            for(unsigned i=0; i<beads.n; i++){
                auto src=beads.get_bead(i);
                unsigned bead_index=bead_hash_to_bead_index.at(src.hash);
                Bead &dst=m_state->beads.at( bead_index );
                for(int d=0; d<3; d++){
                    dst.x[d]=src.x[d];
                    dst.v[d]=src.v[d];
                    dst.f[d]=src.f[d];
                }
                dst.v += dst.f * dt_2;
            }
        }
    }
public:
    void Attach(WorldState *state)
    {
        m_state=state;

        if(!state){
            m_cells.clear();
            return;
        }

        vec3i_t box{(int)state->box[0], (int)state->box[1], (int)state->box[2]};

        const bool use_morton = false;

        unsigned volume;
        if(use_morton){
            volume = m_morton(box[2], box[1], box[0]);
        }else{
            volume = box[0]*box[1]*box[2];
        }

        m_cells.assign(volume, {});

        unsigned offset=0;
        for_each_point_in_box({0,0,0}, box, [&](vec3i_t pos){
            unsigned index= use_morton ? m_morton(pos.x[2], pos.x[1], pos.x[0]) : offset;

            assert(m_cells[index].beads.n==0);
            
            pos_to_index[pos]=index;
            m_cells.at(index).pos=pos;
            m_cells.at(index).is_edge=false;
            ++offset;
        });

        auto rel_nhood=make_relative_nhood_forwards(true);
        rel_nhood.insert(rel_nhood.begin(), vec3i_t(0,0,0)); // Ensure centre is first
        for_each_point_in_box({0,0,0}, box, [&](vec3i_t pos){
            unsigned offset=pos_to_index.at(pos);
            Cell &cell=m_cells[offset];

            auto p=cell.pos;
            unsigned nhood_offset=0;
            for(auto n : make_absolute_nhood(rel_nhood, box, cell.pos)){
                cell.nhood[nhood_offset].reset(&m_cells[pos_to_index[n]].beads);

                for(int d=0; d<3; d++){
                    if(n[d] == box[d]-1 && p[d]==0){
                        cell.nhood[nhood_offset].add_wrap_bits( 2<<(2*d) );
                        cell.is_edge=true;
                    }
                    if(n[d]==0 && p[d] == box[d]-1){
                        cell.nhood[nhood_offset].add_wrap_bits( 1<<(2*d) );
                        cell.is_edge=true;
                    }
                }
                nhood_offset++;
            }

            assert(cell.nhood[0].get_neighbour() == &cell.beads);
        });

        m_maxBeadTypes=m_state->bead_types.size();
        if(m_maxBeadTypes > 8){
            throw std::runtime_error("This engine can't support more than 8 bead types.");
        }else if(m_maxBeadTypes <=4 ){
            m_maxBeadTypes=4;
        }else{
            m_maxBeadTypes=8;
        }

        m_sharedSqrtDissStrength=true;
        // Put them in aligned memory suutable up to AVX512
        m_conStrength.reset((float*)aligned_alloc(64, sizeof(float)*m_maxBeadTypes*m_maxBeadTypes), [](float *p){ free(p); });
        m_sqrtDissStrength.reset((float*)aligned_alloc(64, sizeof(float)*m_maxBeadTypes*m_maxBeadTypes), [](float *p){ free(p); });
        for(unsigned i=0; i<state->bead_types.size(); i++){
            for(unsigned j=0; j<state->bead_types.size(); j++){
                unsigned dst=i*m_maxBeadTypes+j;
                m_conStrength.get()[dst]=state->interactions[i*state->bead_types.size()+j].conservative;
                m_sqrtDissStrength.get()[dst]=sqrtf(state->interactions[i*state->bead_types.size()+j].dissipative);
                if(m_sqrtDissStrength.get()[dst]!=m_sqrtDissStrength.get()[0]){
                    m_sharedSqrtDissStrength=false;
                }
            }   
        }

        m_rng_scale_s32_to_u = pow_half(24*dpd_maths_core_half_step::kT / m_state->dt) * ldexp(1.0, -31);

        m_working.reset(new vector_packed_beads);

        import_beads();
    }

    template<class TConfig>
    void RunConfig(unsigned nSteps)
    {
        for(unsigned i=0; i<nSteps; i++){
            step_impl<TConfig>();
        }
    }

    void Run(unsigned nSteps) override
    {
        if(K < 8 && m_state->bead_types.size() <= 4){
            if(m_sharedSqrtDissStrength){
                RunConfig<Config_MB4_Shared>(nSteps);
            }else{
                RunConfig<Config_MB4>(nSteps);
            }              
        }else{
            if(m_sharedSqrtDissStrength){
                RunConfig<Config_MB8_Shared>(nSteps);
            }else{
                RunConfig<Config_MB8>(nSteps);
            }            
        }
        export_beads();
    }

private:
    template<class TConfig>
    void step_impl()
    {
        m_t_hash = get_t_hash(m_state->t, m_state->seed);

        vec3i_t box{(int)m_state->box[0], (int)m_state->box[1], (int)m_state->box[2]};
        vec3f_t boxf{box};

        float dt=m_state->dt;

        std::vector<uint8_t> outgoing;
        for(Cell &c : m_cells){
            float origin[3];
            c.pos.extract(origin);

            if(c.beads.n==0){
                continue;
            }

            c.beads.newton(origin, dt, outgoing);

            if(!outgoing.empty()){
                for(int i=outgoing.size()-1; i>=0; i--){
                    auto index=outgoing[i];
                    auto b=c.beads.erase(index);
                    for(int d=0; d<3; d++){
                        int dl=(b.x[d] < 0.0f ), dr = (b.x[d]>=box[d]);
                        int delta=(dl-dr)*box[d];
                        b.x[d] += delta;
                        if(b.x[d] == box[d]){
                            b.x[d]=0; // Check for case where it rounds up
                        }
                        assert(0 <= b.x[d] && b.x[d] < box[d]);
                    }
                    auto dst_pos=vec3_floor(b.x);
                    for(int d=0; d<3; d++){
                        assert(0<=dst_pos[d] && dst_pos[d] < box[d]);
                    }
                    assert(pos_to_index.find(dst_pos)!=pos_to_index.end());
                    auto dst_index=pos_to_index.at(dst_pos);
                    auto &dst_cell=m_cells[dst_index];

                    dst_cell.incoming.push_back(b);
                }
                outgoing.clear();
            }
        }

        for(Cell &c : m_cells){
            flush_incoming(c);
        }

        for(Cell &c : m_cells){
            float origin[3];
            c.pos.extract(origin);

            assert(&c.beads==c.nhood[0].get_neighbour());

            if(c.beads.n==0){
                continue;
            }

            filter_and_calc_force<TConfig>(
                m_rng_scale_s32_to_u,
                m_t_hash,
                &boxf.x[0],
                c.is_edge,

                m_conStrength.get(),
                m_sqrtDissStrength.get(),

                &c.nhood[0],
                m_working.get()
            );

        }

        m_state->t++;
    } 

};

    }; // HYW_NAMESPACE

}; // dpd_maths_highway

HWY_AFTER_NAMESPACE();

#endif
