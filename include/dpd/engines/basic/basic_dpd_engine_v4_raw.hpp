#ifndef basic_dpd_engine_v4_raw_hpp
#define basic_dpd_engine_v4_raw_hpp

#include "dpd/engines/basic/basic_dpd_engine.hpp"

#include "dpd/engines/basic/basic_dpd_engine_v4_raw_handlers.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>

/*
   Reducing the number of barriers:

   - No extra barrier for end of update mom (momentum) update. Pos and mom are merged, and
     we pre-distort and post-correct for the extra (start) and missing (end) mom updates.

   - Sharing and force exchanges occur in the same phase.
*/
class BasicDPDEngineV4Raw
    : public BasicDPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=12;

    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;

    using Handlers = BasicDPDEnginev4RawHandlers<BasicDPDEngineV4Raw>;


    using raw_bead_view_t = Handlers::raw_bead_view_t;
    using raw_angle_bond_info_t = Handlers::raw_angle_bond_info_t;
    using raw_bead_resident_t = Handlers::raw_bead_resident_t;
    using raw_cached_bond_t = Handlers::raw_cached_bond_t;
    using raw_force_input_t = Handlers::raw_force_input_t;
    using device_state_t = Handlers::device_state_t;

    
    std::string CanSupport(const WorldState *s) const override
    {

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
            dst.t=m_state->t;
            dst.inv_root_dt=pow_half(24 * dpd_maths_core_half_step::kT / m_state->dt);
            dst.bond_r0=m_bond_r0;
            dst.bond_kappa=m_bond_kappa;
            for(unsigned i=0; i<m_state->bead_types.size(); i++){
                for(unsigned j=0; j<m_state->bead_types.size(); j++){
                    dst.interactions[i*MAX_BEAD_TYPES+j].conservative=m_state->interactions[i*m_state->bead_types.size()+j].conservative;
                    dst.interactions[i*MAX_BEAD_TYPES+j].sqrt_dissipative=sqrt(m_state->interactions[i*m_state->bead_types.size()+j].dissipative);
                }
            }
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

        for(unsigned i=0; i<nSteps; i++){
            step(states, neighbour_map);
            m_state->t += 1;
        }

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

                #ifndef NDEBUG
                src.checksum=calc_checksum(src);
                #endif

                bead_resident_t tmp;
                Handlers::copy_bead_resident(&tmp, &src);
                dst.beads.push_back(tmp);
            }
        }

        export_beads();
    }

    void step(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map)
    {
        if(ForceLogging::logger()){
            step_impl<true>(states, neighbour_map);
        }else{
            step_impl<false>(states, neighbour_map);
        }
    }


    template<bool EnableLogging>
    void step_impl(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map)
    {
        if(EnableLogging && ForceLogging::logger()){
            ForceLogging::logger()->SetTime(m_state->t);
        }

        //////////////////////////////////////////////////////////////
        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &c : states){
            Handlers::on_barrier(c);
        }

        // Do any migration needed
        for(auto &cell : states){
            raw_bead_resident_t buffer;
            uint32_t rts;
            while( (rts = Handlers::calc_rts(cell)) ){
                assert(rts==Handlers::RTS_FLAG_migrate);
                Handlers::on_send_migrate(cell, buffer);
                for(auto ni : neighbour_map[&cell]){
                    Handlers::on_recv_migrate(*ni, buffer);
                }
            }
        }

        ///////////////////////////////////////////////////////////////
        // Shared the beads and calculate DPD and bonded forces
        for(auto &c : states){
            Handlers::on_barrier(c);
        }

        // Calculate all the shared forces, and also do angle forces when able
        bool idle=false;
        while(!idle){
            idle=true;
            for(auto &cell : states){
                uint32_t rts;
                if( (rts = Handlers::calc_rts(cell)) ){
                    if(rts&Handlers::RTS_FLAG_force){
                        raw_force_input_t transfer;
                        assert(cell.force_outgoing.n>0);
                        Handlers::on_send_force(cell, transfer);

                        for(auto ni : neighbour_map[&cell]){
                            Handlers::on_recv_force(*ni, transfer);
                        }
                        idle=false;
                    }else if(rts&Handlers::RTS_FLAG_share){
                        raw_bead_view_t transfer;
                        Handlers::on_send_share(cell, transfer);
                        for(auto ni : neighbour_map[&cell]){
                            Handlers::on_recv_share<EnableLogging>(*ni, transfer);
                        }
                        idle=false;
                    }else{
                        assert(false);
                    }
                }
            }
        }

        if(EnableLogging && ForceLogging::logger()){
            for(auto &c : m_cells){
                for(auto &b : c.beads){
                    double h=b.get_hash_code().hash;
                    double x[3]={b.x[0],b.x[1],b.x[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x_next",3,x);
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_next",3,f);
                }
            }
        }
    }


};

#endif
