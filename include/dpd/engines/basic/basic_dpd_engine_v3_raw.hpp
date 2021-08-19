#ifndef basic_dpd_engine_v3_raw_hpp
#define basic_dpd_engine_v3_raw_hpp

#include "dpd/core/basic_dpd_engine.hpp"

#include "dpd/engines/basic/basic_dpd_engine_v3_raw_handlers.hpp"

#include <iterator>
#include <unordered_map>
#include <cstring>

/*
    Lowering the 
*/
class BasicDPDEngineV3Raw
    : public BasicDPDEngine
{
public:

    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=12;

    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;

    using Handlers = BasicDPDEngineV3RawHandlers<BasicDPDEngineV3Raw>;


    using raw_bead_view_t = Handlers::raw_bead_view_t;
    using raw_angle_bond_info_t = Handlers::raw_angle_bond_info_t;
    using raw_bead_resident_t = Handlers::raw_bead_resident_t;
    using raw_cached_bond_t = Handlers::raw_cached_bond_t;
    using raw_force_input_t = Handlers::raw_force_input_t;
    using bead_bag = Handlers::bead_bag;
    using cached_bond_bag = Handlers::cached_bond_bag;
    using force_input_bag = Handlers::force_input_bag;
    using device_state_t = Handlers::device_state_t;

    
    std::string CanSupport(const WorldState *s) const override
    {
        if(s->bead_types.size()>1){
            double diss=s->interactions.at(1).dissipative;
            for(int i=0; i<s->bead_types.size(); i++){
                for(int j=0; j<s->bead_types.size(); j++){
                    if(s->interactions[i*s->bead_types.size()+j].dissipative!=diss){
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
            dst.phase=Handlers::PreMigrate;

            auto resident=make_bag_wrapper(dst.resident);
            for(const auto &b : src.beads){
                raw_bead_resident_t tmp;
                Handlers::copy_bead_resident(&tmp, &b);
                resident.push_back(tmp);
            }

            for(unsigned neighbour_index : src.neighbours){
                neighbour_map[&dst].push_back( &states[neighbour_index] );
            }
        
        }

        for(unsigned i=0; i<nSteps; i++){
            step(states, neighbour_map);
        }

        //////////////////////////////////////////////
        // Export shadow back out
        for(unsigned i=0; i<m_cells.size(); i++){
            const device_state_t &src=states[i];
            cell_t &dst=m_cells[i];

            auto resident=make_bag_wrapper(src.resident);
            dst.beads.clear();
            for(const auto &src : resident){
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

    void step(std::vector<device_state_t> &states, std::unordered_map<device_state_t*,std::vector<device_state_t*>> &neighbour_map)
    {
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

        // Calculate all the non-angle forces
        for(auto &cell : states){
            raw_bead_view_t transfer;
            uint32_t rts;
            while( (rts = Handlers::calc_rts(cell)) ){
                assert(rts==Handlers::RTS_FLAG_share);
                Handlers::on_send_share(cell, transfer);

                for(auto ni : neighbour_map[&cell]){
                    Handlers::on_recv_share(*ni, transfer);
                }
            }
        }

        ///////////////////////////////////////////////////////////////
        // Shared the beads, now need to do angle forces
        for(auto &cell : states){
            Handlers::on_barrier(cell);
        }

        for(auto &cell : states){
            raw_force_input_t transfer;
            uint32_t rts;
            while( (rts = Handlers::calc_rts(cell)) ){
                assert(rts==Handlers::RTS_FLAG_force);
                assert(cell.force_outgoing.n>0);
                Handlers::on_send_force(cell, transfer);

                for(auto ni : neighbour_map[&cell]){
                    Handlers::on_recv_force(*ni, transfer);
                }
            }
        }

        //////////////////////////////////////////
        // Finally do the velocity updates
        for(auto &cell : states){
            Handlers::on_barrier(cell);
        }
    }


};

#endif
