#ifndef basic_dpd_engine_v2_hpp
#define basic_dpd_engine_v2_hpp

#include "dpd/engines/basic/basic_dpd_engine.hpp"



/*
    This is the same as BasicDPDEngine, but uses a more
    message-like architecture.

    It also uses three barriers, when only 2 are needed.
    The pre-step position and post-step momentum updates could
    be merged, as long as care is taken at the beginning and
    end of a sequence of steps.
*/
class BasicDPDEngineV2
    : public BasicDPDEngine
{
private:

    void on_barrier_pre_migrate(cell_t &c)
    {
        c.cached_bonds.clear();

        unsigned i=0;
        while(i<c.beads.size()){
            auto &b = c.beads[i];
            // Actual position update, clears force
            dpd_maths_core_half_step::update_pos(m_dt, m_box, b);

            // Check it is still in the right cell
            vec3i_t true_loc=floor(b.x);
            if(true_loc != c.location){
                // Movement is fairly unlikely
                c.outgoing.push_back(b);
                if(i+1<c.beads.size()){
                    std::swap(c.beads[i], c.beads.back());
                }
                c.beads.resize(c.beads.size()-1);
            }else{
                ++i;
            }
        }
    }

    void on_send_migrate(cell_t &cell, std::vector<bead_resident_t> &outgoing)
    {
        outgoing.clear();
        std::swap(outgoing, cell.outgoing);
    }

    void on_recv_migrate(cell_t &cell, std::vector<bead_resident_t> &incoming)
    {
        for(const auto &b : incoming){
            if( floor(b.x) == cell.location ){
                cell.beads.push_back(b);
            }
        }
    }

    void on_barrier_pre_share(cell_t &)
    {}

    void on_send_share(cell_t &cell, std::vector<bead_view_t> &outgoing)
    {
        outgoing.clear();
        for(const auto &b : cell.beads){
            outgoing.push_back({ b.id, b.x, b.v });
        }
    }

    template<bool EnableLogging>
    void on_recv_share(cell_t &cell, std::vector<bead_view_t> &incoming)
    {
        // Stay on one neighbour so we can tell if it must be cached
        // We do wrapping on a per bead basis, to simulate receiving each one
        for(const bead_view_t &neighbour_bead : incoming){
            bool cache_neighbour=false;

            vec3f_t neighbour_x = neighbour_bead.x;
            vec3i_t neighbour_cell_pos = floor(neighbour_x);
            for(int d=0; d<3; d++){
                if(cell.location[d]==0 && neighbour_cell_pos[d]==m_box[d]-1){
                    neighbour_x[d] -= m_box[d];
                }else if(cell.location[d]==m_box[d]-1 && neighbour_cell_pos[d]==0){
                    neighbour_x[d] += m_box[d];
                }
            }

            for(auto &bead : cell.beads){
                // This implicitly interacts each bead with itself, which is handled with a
                // distance check in calc_force.
                cache_neighbour |= calc_force<EnableLogging>( bead, neighbour_bead, neighbour_x);
            }

            if(cache_neighbour){
                cell.cached_bonds.push_back({neighbour_bead.id.get_hash_code(), neighbour_x});
            }
        }
    }

    template<bool EnableLogging>
    void on_send_force(cell_t &cell, std::vector<force_input_t> &outgoing)
    {
        for(auto &b : cell.beads){
            update_bead_angle<EnableLogging>(cell.cached_bonds, b, outgoing);
        }
    }

    void on_recv_force(cell_t &cell, std::vector<force_input_t> &outgoing)
    {
        for(auto &b : cell.beads){
            if(!b.id.is_monomer()){
                auto bpo = b.id.get_hash_code();
                for(const auto &f : outgoing){
                    if( bead_hash_equals(bpo , f.target_hash)){
                        b.f += f.f;
                    }
                }
            }
        }
    }

    void on_barrier_post_share(cell_t &cell)
    {
        for(auto &b : cell.beads){
            dpd_maths_core_half_step::update_mom(m_dt, b);
            #ifndef NDEBUG
            b.checksum=calc_checksum(b);
            #endif
        }
    }

    void step() override
    {
        if(ForceLogging::logger()){
            step_impl<true>();
        }else{
            step_impl<false>();
        }
    }

    
    template<bool EnableLogging>
    void step_impl()
    {
        double dt=m_state->dt;
        m_t_hash = get_t_hash(m_state->t, m_state->seed);

        if(EnableLogging && ForceLogging::logger()){
            ForceLogging::bead_hash_to_id() = [&](uint32_t hash) -> uint32_t { return m_state->bead_hash_to_id(hash); };

            ForceLogging::logger()->SetTime(m_state->t);
            ForceLogging::logger()->LogProperty("dt", 1, &m_state->dt);
            double seed_low=m_state->seed &0xFFFFFFFFul;
            double seed_high=m_state->seed>>32;
            ForceLogging::logger()->LogProperty("seed_lo", 1, &seed_low);
            ForceLogging::logger()->LogProperty("seed_high", 1, &seed_high);
            double t_hash_low=m_t_hash&0xFFFFFFFFul;
            double t_hash_high=m_t_hash>>32;
            ForceLogging::logger()->LogProperty("t_hash_lo", 1, &t_hash_low);
            ForceLogging::logger()->LogProperty("t_hash_high", 1, &t_hash_high);
            for(auto &c : m_cells){
                for(auto &b : c.beads){
                    double h=b.get_hash_code();
                    auto bead_id=m_state->polymers.at(b.id.get_polymer_id()).bead_ids.at(b.id.get_polymer_offset());
                    ForceLogging::logger()->LogBeadProperty(bead_id, "b_hash", 1, &h);
                    double x[3]={b.x[0],b.x[1],b.x[2]};
                    ForceLogging::logger()->LogBeadProperty(bead_id,"x",3,x);
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(bead_id,"f",3,f);
                }
            }
        }

        //////////////////////////////////////////////////////////////
        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &c : m_cells){
            on_barrier_pre_migrate(c);
        }

        // Do any migration needed
        for(auto &c : m_cells){
            std::vector<bead_resident_t> buffer;
            on_send_migrate(c, buffer);
            for(auto ni : c.neighbours){
                on_recv_migrate(m_cells[ni], buffer);
            }
        }

        ///////////////////////////////////////////////////////////////
        // Shared the beads and calculate DPD and bonded forces
        for(auto &c : m_cells){
            on_barrier_pre_share(c);
        }

        // Calculate all the non-angle forces
        for(cell_t &cell : m_cells){
            std::vector<bead_view_t> transfer;
            on_send_share(cell, transfer);
            for(unsigned neighbour_index : cell.neighbours){
                on_recv_share<EnableLogging>(m_cells[neighbour_index], transfer);
            }
        }

        // At some point we have received and cached all angle bond updates.
        // An extra mechanism is needed to know when this can be done
        // This must be done pre hardware idle, as any forces need to be send out
        for(cell_t &cell : m_cells){
            std::vector<force_input_t> transfer;
            on_send_force<EnableLogging>(cell, transfer);
            for(unsigned neighbour_index : cell.neighbours){
                on_recv_force(m_cells[neighbour_index], transfer);
            }
        }

        //////////////////////////////////////////
        // Finally do the velocity updates
        for(cell_t &cell : m_cells){
            on_barrier_post_share(cell);
        }

        m_state->t += 1;

        if(EnableLogging && ForceLogging::logger()){
            for(auto &c : m_cells){
                for(auto &b : c.beads){
                    double h=b.get_hash_code();
                    auto bead_id=m_state->polymers.at(b.id.get_polymer_id()).bead_ids.at(b.id.get_polymer_offset());
                    double x[3]={b.x[0],b.x[1],b.x[2]};
                    ForceLogging::logger()->LogBeadProperty(bead_id,"x_next",3,x);
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(bead_id,"f_next",3,f);
                }
            }
            ForceLogging::bead_hash_to_id() = {};
        }
    }


};

#endif
