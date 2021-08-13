#ifndef basic_dpd_engine_v3_hpp
#define basic_dpd_engine_v3_hpp

#include "basic_dpd_engine.hpp"
#include "storage/bag_wrapper.hpp"

#include <iterator>


/*
    This is even more message-like than BasicDPDEngineV3

    For convenience this has risen to four barriers, to ensure that all
    cached bonds are seen before forces are recieved.
*/
class BasicDPDEngineV3
    : public BasicDPDEngine
{
private:
    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    struct bead_bag
    {
        bead_resident_t elements[MAX_BEADS_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };

    struct cached_bond_bag
    {
        cached_bond_t elements[MAX_CACHED_BONDS_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };

    struct force_input_bag
    {
        force_input_t elements[MAX_OUTGOING_FORCES_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };

    struct device_state_t
    {
        int32_t box[3];
        float dt;

        int32_t location[3];

        uint32_t share_todo;
        bead_bag resident;
        bead_bag migrate_outgoing;
        force_input_bag force_outgoing;
        cached_bond_bag cached_bonds;
    };


    void on_barrier_pre_migrate(device_state_t &c)
    {
        auto cached_bonds=make_bag_wrapper(c.cached_bonds);
        auto resident=make_bag_wrapper(c.resident);
        auto migrate_outgoing=make_bag_wrapper(c.migrate_outgoing);

        cached_bonds.clear();

        assert(migrate_outgoing.empty());

        float dt=c.dt;
        vec3i_t box{c.box};
        vec3i_t location{c.location};
        
        int i=resident.size();
        while(0 < i){
            --i;
            auto &b = resident[i];
            // Actual position update, clears force
            dpd_maths_core_half_step::update_pos(dt, box, b);
        
            // Check it is still in the right cell
            vec3i_t true_loc=floor(b.x);
            if(true_loc != c.location){
                // Movement is fairly unlikely
                migrate_outgoing.push_back(b);

                std::swap(resident.back(), b); // A waste if this is the last bead, but makes code simpler and smaller
                resident.pop_back();
            }
        }
    }

    bool is_rts_migrate(const device_state_t &c)
    { return !make_bag_wrapper(c.migrate_outgoing).empty(); }

    void on_send_migrate(device_state_t &c, bead_resident_t &outgoing)
    {
        auto migrate_outgoing=make_bag_wrapper(c.migrate_outgoing);
        assert(!migrate_outgoing.empty());
        
        outgoing=migrate_outgoing.back();
        migrate_outgoing.pop_back();
    }

    void on_recv_migrate(device_state_t &c, bead_resident_t &incoming)
    {
        auto resident=make_bag_wrapper(c.resident);

        vec3i_t incoming_loc=floor(incoming.x);
        if(incoming_loc == vec3i_t(c.location)){
            resident.push_back(incoming);
        }
    }

    void on_barrier_pre_share(device_state_t &c)
    {
        auto resident=make_bag_wrapper(c.resident);
        c.share_todo = resident.size();
    }

    bool is_rts_share(const device_state_t &c)
    { return c.share_todo>0; }

    void on_send_share(device_state_t &cell, bead_view_t &outgoing)
    {
        auto resident=make_bag_wrapper(cell.resident);

        assert(cell.share_todo>0);
        --cell.share_todo;
        const auto &b = resident[cell.share_todo];
        outgoing.id = b.id;
        outgoing.x = b.x;
        outgoing.v = b.v;
    }

    void on_recv_share(device_state_t &cell, const bead_view_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);

        vec3i_t box{cell.box};
        
        vec3f_t neighbour_x{incoming.x};
        vec3i_t neighbour_cell_pos = floor(neighbour_x);
        for(int d=0; d<3; d++){
            if(cell.location[d]==0 && neighbour_cell_pos[d]==m_box[d]-1){
                neighbour_x[d] -= m_box[d];
            }else if(cell.location[d]==m_box[d]-1 && neighbour_cell_pos[d]==0){
                neighbour_x[d] += m_box[d];
            }
        }

        bool cache_neighbour=false;


        for(auto &bead : resident){
            // This implicitly interacts each bead with itself, which is handled with a
            // distance check in calc_force.
            cache_neighbour |= calc_force( bead, incoming, neighbour_x);
        }

        if(cache_neighbour){
            cached_bonds.push_back({incoming.id.get_hash_code(), neighbour_x});
        }
    }

    void on_barrier_pre_force(device_state_t &cell)
    {
        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);
        assert(force_outgoing.empty());

        for(auto &b : resident){
            update_bead_angle(cached_bonds, b, force_outgoing);
        }

        cached_bonds.clear();
    }

    bool is_rts_force(const device_state_t &c)
    { return !make_bag_wrapper(c.force_outgoing).empty(); }

    void on_send_force(device_state_t &cell,force_input_t &outgoing)
    {
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

        assert(!force_outgoing.empty());
        outgoing=force_outgoing.back();
        force_outgoing.pop_back();
    }

    void on_recv_force(device_state_t &cell, const force_input_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        for(auto &b : resident){
            if(b.id.get_hash_code() == incoming.target_hash){
                b.f += incoming.f;
            }
        }
    }

    void on_barrier_post_share(device_state_t &cell)
    {
        auto resident=make_bag_wrapper(cell.resident);

        for(auto &b : resident){
            dpd_maths_core_half_step::update_mom(cell.dt, b);
        }
    }
    

    void step() override
    {
        double dt=m_state->dt;
        m_t_hash = next_t_hash(m_state->seed);

        // Create shadow device states
        std::vector<device_state_t> states;
        std::unordered_map<device_state_t*,std::vector<device_state_t*>> neighbour_map;
        states.resize(m_cells.size());
        for(unsigned i=0; i<m_cells.size(); i++){
            device_state_t &dst=states[i];
            const cell_t &src=m_cells[i];

            dst.dt=dt;
            m_box.extract(dst.box);
            src.location.extract(dst.location);

            auto resident=make_bag_wrapper(dst.resident);
            for(const auto &b : src.beads){
                resident.push_back(b);
            }

            for(unsigned neighbour_index : src.neighbours){
                neighbour_map[&dst].push_back( &states[neighbour_index] );
            }
        }


        //////////////////////////////////////////////////////////////
        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &c : states){
            on_barrier_pre_migrate(c);
        }

        // Do any migration needed
        for(auto &c : states){
            bead_resident_t buffer;
            while(is_rts_migrate(c)){
                on_send_migrate(c, buffer);
                for(auto ni : neighbour_map[&c]){
                    on_recv_migrate(*ni, buffer);
                }
            }
        }

        ///////////////////////////////////////////////////////////////
        // Shared the beads and calculate DPD and bonded forces
        for(auto &c : states){
            on_barrier_pre_share(c);
        }

        // Calculate all the non-angle forces
        for(auto &cell : states){
            bead_view_t transfer;
            while(is_rts_share(cell)){
                on_send_share(cell, transfer);

                for(auto ni : neighbour_map[&cell]){
                    on_recv_share(*ni, transfer);
                }
            }
        }

        ///////////////////////////////////////////////////////////////
        // Shared the beads, now need to do angle forces
        for(auto &cell : states){
            on_barrier_pre_force(cell);
        }

        for(auto &cell : states){
            force_input_t transfer;
            while(is_rts_force(cell)){
                on_send_force(cell, transfer);

                for(auto ni : neighbour_map[&cell]){
                    on_recv_force(*ni, transfer);
                }
            }
        }

        //////////////////////////////////////////
        // Finally do the velocity updates
        for(auto &cell : states){
            on_barrier_post_share(cell);
        }

        //////////////////////////////////////////////
        // Export back out

        for(unsigned i=0; i<m_cells.size(); i++){
            const device_state_t &src=states[i];
            cell_t &dst=m_cells[i];

            auto resident=make_bag_wrapper(src.resident);
            dst.beads.assign(resident.begin(), resident.end());
        }

        m_state->t += m_state->dt;
    }


};

#endif
