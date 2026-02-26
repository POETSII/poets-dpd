#ifndef dpd_engine_highway_non_conflict_hpp
#define dpd_engine_highway_non_conflict_hpp

#include "dpd/engines/highway/dpd_engine_highway_single_thread.hpp"

#include "dpd/core/non_conflict_grid.hpp"

HWY_BEFORE_NAMESPACE();

namespace dpd_maths_highway
{
    namespace HWY_NAMESPACE
    {
        namespace hn = hwy::HWY_NAMESPACE;

template<bool TUseHash=false>
class DPDEngineHighwayNonConflict
    : public DPDEngingHighwaySingleThread<TUseHash>
{
private:

public:
    void Attach(WorldState *state)
    {
        DPDEngingHighwaySingleThread<TUseHash>::Attach(state);
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
        if(m_state->bead_types.size() <= 4){
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
                    auto dst_pos=vec3_floor(b.x);
                    for(int d=0; d<3; d++){
                        int dl=(dst_pos[d]<0), dr=(dst_pos[d]>=box[d]);
                        int delta=(dl - dr) * box[d];
                        dst_pos[d] +=delta;
                        b.x[d] += (float) delta;
                        assert(0<=dst_pos[d] && dst_pos[d] < box[d]);
                    }
                    assert(pos_to_index.find(dst_pos)!=pos_to_index.end());
                    auto dst_index=pos_to_index.at(dst_pos);
                    m_cells[dst_index].incoming.push_back(b);
                }
                outgoing.clear();
            }
        }

        // TODO!! : Don't do this in a seperate loop
        for(Cell &c : m_cells){
            if(!c.incoming.empty()){
                // Should enter this infrequently
                flush_incoming(c);
            }
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

                m_conStrength.get(),
                m_sqrtDissStrength.get(),

                c.nhood.size(),
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
