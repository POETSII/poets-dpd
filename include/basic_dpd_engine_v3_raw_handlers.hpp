#ifndef basic_dpd_engine_v3_raw_handlers_hpp
#define basic_dpd_engine_v3_raw_handlers_hpp

#include "storage/bag_wrapper.hpp"
#include "dpd_maths_core_half_step_raw.hpp"

#include <iterator>
#include <cstring>

struct BasicDPDEngineV3RawConfig
{
    static constexpr size_t MAX_BONDS_PER_BEAD = 3;
    static constexpr size_t MAX_BEADS_PER_CELL = 8;
    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=3;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=8;
};

/*
    Lowering the 
*/
template<class TConfig>
struct BasicDPDEngineV3RawHandlers
{
    static constexpr int MAX_BONDS_PER_BEAD = TConfig::MAX_BONDS_PER_BEAD;
    static constexpr int MAX_BEADS_PER_CELL = TConfig::MAX_BEADS_PER_CELL;
    static constexpr int MAX_CACHED_BONDS_PER_CELL = TConfig::MAX_CACHED_BONDS_PER_CELL;
    static constexpr int MAX_OUTGOING_FORCES_PER_CELL = TConfig::MAX_OUTGOING_FORCES_PER_CELL;
    static constexpr int MAX_ANGLE_BONDS_PER_BEAD = TConfig::MAX_ANGLE_BONDS_PER_BEAD;

    static constexpr int MAX_BEAD_TYPES = TConfig::MAX_BEAD_TYPES;

    enum OutputFlags
    {
        RTS_FLAG_migrate = 0x1,
        RTS_FLAG_share = 0x2,
        RTS_FLAG_force = 0x4
    };

    enum Phase
    {
        PreMigrate,
        Migrating,
        Sharing,
        Forcing
    };

    static uint32_t get_bead_type(uint32_t bead_id)
    { return bead_id & 0xF; }

    static bool is_monomer(uint32_t bead_id) 
    { return bead_id>>31; }

    static uint32_t get_polymer_id(uint32_t bead_id)
    { return (bead_id>>4) & (is_monomer(bead_id) ? 0x07FF'FFFFul : 0x000F'FFFFul); }

    static uint32_t get_polymer_offset(uint32_t bead_id)
    { return is_monomer(bead_id) ? 0 : (bead_id>>24)&0x7F; }

    static uint32_t make_hash_from_offset(uint32_t bead_id, unsigned offset)
    {
        assert(!is_monomer(bead_id));
        return ((bead_id>>4)&0x000F'FFFFul) | (offset<<20);
    }

    static uint32_t get_hash_code(uint32_t bead_id)
    { return bead_id>>4; }


    struct raw_bead_view_t
    {
        uint32_t id;
        float x[3];
        float v[3];
    };
    static_assert(sizeof(raw_bead_view_t)==28);

    static uint32_t get_hash_code(const raw_bead_view_t &bead)
    { return get_hash_code(bead.id); }

    static uint32_t get_bead_type(const raw_bead_view_t &bead)
    { return get_bead_type(bead.id); }

    struct raw_angle_bond_info_t
    {
        uint8_t partner_head;
        uint8_t partner_tail;
        uint8_t kappa;   // Kappa is just an integer. Typically this is quite small, e.g. 5 or 15. Should be less than 255
        uint8_t _pad_;
    };

    struct raw_bead_resident_t
    {
        // raw_bead_view_t
        uint32_t id;
        float x[3];
        float v[3];

        // Extras
        float f[3];

        uint8_t bond_partners[MAX_BONDS_PER_BEAD+1]; // -1 means no bond
        static_assert(sizeof(bond_partners)==4); // Keep aligned

        raw_angle_bond_info_t angle_bonds[MAX_ANGLE_BONDS_PER_BEAD];
    };
    static_assert(sizeof(raw_bead_resident_t)==48);

    template<class A, class B>
    static void copy_bead_view(A *dst, const B *src)
    {
        static_assert(offsetof(A,id)==0);
        static_assert(offsetof(A,id)==offsetof(B,id));
        static_assert(offsetof(A,x)==sizeof(A::id));
        static_assert(offsetof(A,x)==offsetof(B,x));
        static_assert(offsetof(A,v)==offsetof(A,x)+sizeof(A::x));
        static_assert(offsetof(A,v)==offsetof(B,v));
        static_assert( (offsetof(A,v)+sizeof(A::v)) % 32 == 4);
        memcpy32((uint32_t*)dst, (const uint32_t*)src, (offsetof(A,v)+sizeof(A::v))/4);
    }


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winvalid-offsetof"
    // TODO : g++ correctly complains, as BasicDPDEngine::bead_resident_t is not strictly POD
    template<class A, class B>
    static void copy_bead_resident(A *dst, const B *src)
    {
        static_assert(sizeof(A)==sizeof(B));
        static_assert(offsetof(A,id)==offsetof(B,id));
        static_assert(offsetof(A,x)==offsetof(B,x));
        static_assert(offsetof(A,v)==offsetof(B,v));
        static_assert(offsetof(A,f)==offsetof(B,f));
        static_assert(offsetof(A,bond_partners)==offsetof(B,bond_partners));
        static_assert(offsetof(A,angle_bonds)==offsetof(B,angle_bonds));
        memcpy32((uint32_t*)dst, (const uint32_t*)src, sizeof(A)/4);
    }
#pragma GCC diagnostic pop

    struct raw_cached_bond_t
    {
        uint32_t bead_hash;
        float x[3];
    };

    struct raw_force_input_t
    {
        uint32_t target_hash;
        float f[3];
    };


    struct bead_bag
    {
        raw_bead_resident_t elements[MAX_BEADS_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };

    struct cached_bond_bag
    {
        raw_cached_bond_t elements[MAX_CACHED_BONDS_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };

    struct force_input_bag
    {
        raw_force_input_t elements[MAX_OUTGOING_FORCES_PER_CELL];
        uint16_t n=0;
        uint16_t lost=0;
    };



    struct device_state_t
    {

        int32_t box[3];
        float dt=nanf("");
        float inv_root_dt=nanf("");
        float bond_r0 = nanf("");
        float bond_kappa = nanf("");
        float conservative[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
        float sqrt_dissipative;
        uint64_t t_hash;
        uint64_t t_seed;

        int32_t location[3];

        Phase phase;
        uint32_t rts;
        uint32_t steps_todo;

        uint32_t share_todo;
        bead_bag resident;
        bead_bag migrate_outgoing;
        force_input_bag force_outgoing;
        cached_bond_bag cached_bonds;
    };

    static uint32_t calc_rts(const device_state_t &cell)
    {
        return cell.rts;
    }

    static void on_barrier(device_state_t &cell)
    {
        switch(cell.phase){
            case PreMigrate: on_barrier_pre_migrate(cell); break;
            case Migrating: on_barrier_pre_share(cell); break;
            case Sharing: on_barrier_pre_force(cell); break;
            case Forcing: on_barrier_post_force(cell); break;
            default: assert(false); break;
        }
    }


    static void on_barrier_pre_migrate(device_state_t &cell)
    {
        assert(cell.phase==PreMigrate);

        auto resident=make_bag_wrapper(cell.resident);
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);

        assert(make_bag_wrapper(cell.cached_bonds).empty());
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        assert(migrate_outgoing.empty());

        int i=resident.size();
        while(0 < i){
            --i;
            auto &b = resident[i];
            // Actual position update, clears force
            dpd_maths_core_half_step_raw::update_pos(cell.dt, cell.box, b);
        
            // Check it is still in the right cell
            int true_loc[3];
            vec3_floor(true_loc, b.x);
            if(true_loc != cell.location){
                // Movement is fairly unlikely
                migrate_outgoing.push_back(b);

                std::swap(resident.back(), b); // A waste if this is the last bead, but makes code simpler and smaller
                resident.pop_back();
            }
        }

        cell.t_hash = next_t_hash(cell.t_seed);

        cell.phase=Migrating;
        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;
    }


    static void on_send_migrate(device_state_t &cell, raw_bead_resident_t &outgoing)
    {
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);
        assert(!migrate_outgoing.empty());
        
        outgoing=migrate_outgoing.back();
        migrate_outgoing.pop_back();

        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;
    }

    static void on_recv_migrate(device_state_t &cell, raw_bead_resident_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        int incoming_loc[3];
        vec3_floor(incoming_loc, incoming.x);
        if(vec3_equal(incoming_loc , cell.location)){
            resident.push_back(incoming);
        }
    }

    static void on_barrier_pre_share(device_state_t &cell)
    {
        assert(cell.phase==Migrating);

        auto resident=make_bag_wrapper(cell.resident);
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        cell.share_todo = resident.size();

        cell.phase=Sharing;
        cell.rts=cell.share_todo==0 ? 0 : RTS_FLAG_share;
    }


    static void on_send_share(device_state_t &cell, raw_bead_view_t &outgoing)
    {
        auto resident=make_bag_wrapper(cell.resident);

        assert(cell.share_todo>0);
        --cell.share_todo;
        const auto &b = resident[cell.share_todo];
        copy_bead_view( &outgoing, &b );
        outgoing.id = b.id;

        cell.rts=cell.share_todo==0 ? 0 : RTS_FLAG_share;
    }

    static void on_recv_share(device_state_t &cell, const raw_bead_view_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        
        float neighbour_x[3];
        vec3_copy(neighbour_x, incoming.x);
        int neighbour_cell_pos[3];
        vec3_floor(neighbour_cell_pos, neighbour_x);
        for(int d=0; d<3; d++){
            if(cell.location[d]==0 && neighbour_cell_pos[d]==cell.box[d]-1){
                neighbour_x[d] -= cell.box[d];
            }else if(cell.location[d]==cell.box[d]-1 && neighbour_cell_pos[d]==0){
                neighbour_x[d] += cell.box[d];
            }
        }

        bool cache_neighbour=false;

        for(auto &bead : resident){
            // This implicitly interacts each bead with itself, which is handled with a
            // distance check in calc_force.
            float dx[3];
            vec3_sub(dx, bead.x, neighbour_x);
            float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            if(dr_sqr >=1 || dr_sqr < 1e-5){ // The min threshold avoid large forces, and also skips self-interaction
                continue;
            }
            float dr=pow_half(dr_sqr);

            bool is_bonded=false;
            float kappa=0.0f;
            float r0=cell.bond_r0;
            if(!(is_monomer(incoming.id) && is_monomer(bead.id) && (get_polymer_id(bead.id) == get_polymer_id(incoming.id)) ) ){
                auto other_polymer_offset=get_polymer_offset(incoming.id);
                for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
                    if(other_polymer_offset==bead.bond_partners[i]){
                        is_bonded=true;
                        kappa=cell.bond_kappa;
                        break;
                    }
                }
            }

            float conStrength=cell.conservative[MAX_BEAD_TYPES*get_bead_type(bead.id)+get_bead_type(incoming.id)];

            float f[3];
            dpd_maths_core_half_step_raw::calc_force<float,float[3],float[3]>(
                cell.inv_root_dt,
                cell.t_hash,
                dx, dr,
                kappa, r0, 
                conStrength,
                cell.sqrt_dissipative,
                get_hash_code(bead.id), get_hash_code(incoming.id),
                bead.v, incoming.v,
                f
            );

            vec3_add(bead.f, f);
            cache_neighbour |= is_bonded;
        }

        if(cache_neighbour){
            raw_cached_bond_t tmp;
            tmp.bead_hash=get_hash_code(incoming.id);
            vec3_copy(tmp.x, neighbour_x);
            cached_bonds.push_back(tmp);
        }
    }

    static void on_barrier_pre_force(device_state_t &cell)
    {
        assert(cell.phase==Sharing);

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);
        assert(force_outgoing.empty());

        auto find_cached_pos=[&](uint32_t target_hash) -> const raw_cached_bond_t *
        {
            for(auto &cb : cached_bonds){
                if(cb.bead_hash==target_hash){
                    return &cb;
                }
            }
            assert(0); // Should be impossible to get here unless bond has snapped
            return 0;
        };

        for(auto &bead : resident){
            for(int i=0; i<MAX_ANGLE_BONDS_PER_BEAD; i++){
                if(bead.angle_bonds[i].partner_head==0xFF){
                    break;
                }
                auto head_hash=make_hash_from_offset(bead.id, bead.angle_bonds[i].partner_head);
                auto tail_hash=make_hash_from_offset(bead.id, bead.angle_bonds[i].partner_tail);
                const auto *head=find_cached_pos(head_hash);
                const auto *tail=find_cached_pos(tail_hash);
                
                // The cache copies should already have wrapping applied
                float first[3], second[3];
                vec3_sub(first, bead.x, head->x);
                vec3_sub(second, tail->x, bead.x);

                float FirstLength= vec3_l2_norm(first);
                float SecondLength  = vec3_l2_norm(second);
                assert(FirstLength < 1);
                assert(SecondLength < 1);

                float headForce[3], middleForce[3], tailForce[3];

                dpd_maths_core_half_step_raw::calc_angle_force<true,float,float[3],float[3]>(
                    (float)bead.angle_bonds[i].kappa, 0.0f, 0.0f,
                    first, FirstLength,
                    second, SecondLength,
                    headForce, middleForce, tailForce
                );

                vec3_add(bead.f, middleForce);

                raw_force_input_t tmp1;
                tmp1.target_hash=head_hash;
                vec3_copy(tmp1.f, headForce);
                force_outgoing.push_back( tmp1 );

                raw_force_input_t tmp2;
                tmp2.target_hash=tail_hash;
                vec3_copy(tmp2.f, tailForce);
                force_outgoing.push_back( tmp2 );
            }
        }

        cached_bonds.clear();

        cell.phase=Forcing;
        cell.rts=force_outgoing.empty() ? 0 : RTS_FLAG_force;
    }


    static void on_send_force(device_state_t &cell,raw_force_input_t &outgoing)
    {
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

        assert(!force_outgoing.empty());
        outgoing=force_outgoing.back();
        force_outgoing.pop_back();

        cell.rts=force_outgoing.empty() ? 0 : RTS_FLAG_force;
    }

    static void on_recv_force(device_state_t &cell, const raw_force_input_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        for(auto &b : resident){
            if( get_hash_code(b.id) == incoming.target_hash){
                vec3_add(b.f, incoming.f);
            }
        }
    }

    static void on_barrier_post_force(device_state_t &cell)
    {
        assert(cell.phase==Forcing);

        auto resident=make_bag_wrapper(cell.resident);
        assert(make_bag_wrapper(cell.cached_bonds).empty());
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        for(auto &b : resident){
            dpd_maths_core_half_step_raw::update_mom(cell.dt, b);
        }

        cell.phase=PreMigrate;
    }

};

#endif
