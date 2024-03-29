#ifndef basic_dpd_engine_v4_raw_handlers_hpp
#define basic_dpd_engine_v4_raw_handlers_hpp

#include "dpd/core/dpd_state.hpp"
#include "dpd/storage/bag_wrapper.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"

#include <iterator>
#include <cstring>

struct BasicDPDEnginev4RawConfig
{
    static constexpr size_t MAX_BONDS_PER_BEAD = 4;
    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=8;
};


template<class TConfig = BasicDPDEnginev4RawConfig>
struct BasicDPDEnginev4RawHandlers
{
    static constexpr int MAX_BONDS_PER_BEAD = TConfig::MAX_BONDS_PER_BEAD;
    static constexpr int MAX_BEADS_PER_CELL = TConfig::MAX_BEADS_PER_CELL;
    static constexpr int MAX_CACHED_BONDS_PER_CELL = TConfig::MAX_CACHED_BONDS_PER_CELL;
    static constexpr int MAX_OUTGOING_FORCES_PER_CELL = TConfig::MAX_OUTGOING_FORCES_PER_CELL;
    static constexpr int MAX_ANGLE_BONDS_PER_BEAD = TConfig::MAX_ANGLE_BONDS_PER_BEAD;

    static constexpr int MAX_BEAD_TYPES = TConfig::MAX_BEAD_TYPES;

    static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);

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
        SharingAndForcing
    };

    static uint32_t get_bead_type(uint32_t bead_id)
    { return BeadHash{bead_id}.get_bead_type(); }

    static bool is_monomer(uint32_t bead_id) 
    { return BeadHash{bead_id}.is_monomer(); }

    static uint32_t get_polymer_id(uint32_t bead_id)
    { return BeadHash{bead_id}.get_polymer_id(); }

    static uint32_t get_polymer_offset(uint32_t bead_id)
    { return BeadHash{bead_id}.get_polymer_offset(); }

    static BeadHash make_hash_from_offset(uint32_t bead_id, unsigned offset)
    {
        return BeadHash{bead_id}.make_reduced_hash_from_polymer_offset( offset);
    }

    static BeadHash get_hash_code(uint32_t bead_id)
    { return BeadHash{bead_id}; }


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

    

    struct raw_bead_resident_t
    {
        // raw_bead_view_t
        uint32_t id;
        float x[3];
        float v[3];

        // Extras
        float f[3];

        uint8_t bond_partners[MAX_BONDS_PER_BEAD]; // -1 means no bond
        static_assert(sizeof(bond_partners)==4); // Keep aligned

        struct raw_angle_bond_info_t
        {
            uint8_t partner_head;
            uint8_t partner_tail;
            uint8_t kappa;   // Kappa is just an integer. Typically this is quite small, e.g. 5 or 15. Should be less than 255
            uint8_t _pad_;
        }angle_bonds[MAX_ANGLE_BONDS_PER_BEAD];

        uint32_t t;
        uint32_t checksum;
    };
    static_assert(sizeof(raw_bead_resident_t)==56);

    using raw_angle_bond_info_t = decltype(raw_bead_resident_t::angle_bonds[0]);

    template<class A, class B>
    static void copy_bead_view(A *dst, const B *src)
    {
        static_assert(offsetof(A,id)==0);
        static_assert(offsetof(A,id)==offsetof(B,id));
        static_assert(offsetof(A,x)==sizeof(A::id));
        static_assert(offsetof(A,x)==offsetof(B,x));
        static_assert(offsetof(A,v)==offsetof(A,x)+sizeof(A::x));
        static_assert(offsetof(A,v)==offsetof(B,v));
        memcpy(dst, src, offsetof(A,v)+sizeof(A::v));
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
        memcpy((char*)dst, (char*)src, sizeof(A));
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

    template<class A, class B>
    static void copy_force_input(A *dst, const B *src)
    {
        static_assert(sizeof(A)==sizeof(B));
        static_assert(offsetof(A,target_hash)==offsetof(B,target_hash));
        static_assert( std::is_same<decltype(A::f),decltype(B::f)>::value );
        memcpy(dst, src, sizeof(B));
    }

    struct device_state_t
    {

        int32_t box[3];
        float dt=nanf("");
        float inv_root_dt=nanf("");
        float bond_r0 = nanf("");
        float bond_kappa = nanf("");
        struct {
            float conservative;
            float sqrt_dissipative;
        }interactions[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
        uint32_t t;
        uint64_t t_hash;
        uint64_t t_seed;

        int32_t location[3];

        Phase phase;
        uint32_t rts;
        uint32_t steps_todo;

        uint32_t share_todo;
        struct{
            raw_bead_resident_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }resident;
        struct{
            raw_bead_resident_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }migrate_outgoing;
        struct force_input_bag
        {
            struct raw_force_input_t
            {
                uint32_t target_hash;
                float f[3];
            } elements[MAX_OUTGOING_FORCES_PER_CELL];
            uint16_t n;
            uint16_t lost;
        } force_outgoing;
        struct cached_bond_bag
        {
            raw_cached_bond_t elements[MAX_CACHED_BONDS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        } cached_bonds;

        // Should probably be co-located with the beads for caching
        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
        uint8_t cached_bond_indices[MAX_BEADS_PER_CELL]; // 0xff if not seen, <0xfe is the first bond that was seen, 0xfe means both beads seen
    
        #ifndef NDEBUG
        uint8_t debug_hookean_bonds_seen[MAX_BEADS_PER_CELL];
        #endif
    };

    static uint32_t calc_rts(const device_state_t &cell)
    {
        return cell.rts;
    }

    template<class TState>
    static void on_barrier(TState &cell)
    {
        switch(cell.phase){
            case PreMigrate:
            case SharingAndForcing: on_barrier_pre_migrate(cell); break;
            case Migrating: on_barrier_pre_share(cell); break;
            default: assert(false); break;
        }
    }


    static void on_barrier_pre_migrate(device_state_t &cell)
    {
        assert(cell.phase==PreMigrate || cell.phase==SharingAndForcing);

        auto resident=make_bag_wrapper(cell.resident);
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);

        make_bag_wrapper(cell.cached_bonds).clear();

        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);
        #ifndef NDEBUG
        if(cell.phase==SharingAndForcing){
            for(unsigned i=0; i<resident.size(); i++){
                const auto &b = resident[i];
                if(b.angle_bonds[0].partner_head!=0xFF){
                    assert(cell.cached_bond_indices[i]==0xFE); // Check angle was processed
                }
                int nhookean=0;
                for(int j=0;j<MAX_BONDS_PER_BEAD;j++){
                    if(b.bond_partners[j]==0xFF){
                        break;
                    }
                    nhookean++;
                }
                assert(cell.debug_hookean_bonds_seen[i]==nhookean); // Check all bonds processed
            }
        }
        #endif

        assert(migrate_outgoing.empty());

        int i=resident.size();
        while(0 < i){
            --i;
            auto &b = resident[i];

            // We always apply the mom step here, rather than doing it in 
            // a seperate step at the end. This means that the input needs
            // to have been "reverse" update_mom'd, and the output must be
            // mom'd again.
            dpd_maths_core_half_step_raw::update_mom(cell.dt, b);

            // Actual position update, clears force
            dpd_maths_core_half_step_raw::update_pos(cell.dt, cell.box, b);

            #ifndef NDEBUG
            b.checksum=calc_checksum(b);
            #endif
        
            // Check it is still in the right cell
            int32_t true_loc[3];
            vec3_floor_nn(true_loc, b.x);
            if( !vec3_equal(true_loc, cell.location) ){
                // Movement is fairly unlikely
                migrate_outgoing.push_back(b);

                std::swap(resident.back(), b); // A waste if this is the last bead, but makes code simpler and smaller
                resident.pop_back();
            }
        }

        cell.t_hash = get_t_hash(cell.t, cell.t_seed);
        cell.t += 1;

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

        int32_t incoming_loc[3];
        vec3_floor_nn(incoming_loc, incoming.x);
        if(vec3_equal(incoming_loc , cell.location)){
            resident.push_back(incoming);
        }else{
            for(int i=0; i<3; i++){
                if(incoming_loc[i] < 0){
                    printf("<m\n");
                }
                if(incoming_loc[i] >= cell.box[i]){
                    printf(">m\n");
                }
            }
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

        for(unsigned i=0; i<resident.size(); i++){ 
            assert( vec3f_t{resident[i].f}.l2_norm() == 0);
            cell.cached_bond_indices[i]=0xff;
            #ifndef NDEBUG
            cell.debug_hookean_bonds_seen[i]=0;
            #endif
        }

        cell.phase=SharingAndForcing;
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

        if(cell.share_todo==0){
            cell.rts &= ~RTS_FLAG_share;
        }
    }

    template<bool EnableLogging>
    static void on_recv_share(device_state_t &cell, const raw_bead_view_t &incoming)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<"), p="<<&cell<<", nres="<<cell.resident.n<<", other="<<get_hash_code(incoming.id)<<"\n";

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);
        
        float neighbour_x[3];
        vec3_copy(neighbour_x, incoming.x);
        int32_t neighbour_cell_pos[3];
        vec3_floor_nn(neighbour_cell_pos, neighbour_x);
        for(int d=0; d<3; d++){
            if(cell.location[d]==0 && neighbour_cell_pos[d]==cell.box[d]-1){
                neighbour_x[d] -= cell.box[d];
            }else if(cell.location[d]==cell.box[d]-1 && neighbour_cell_pos[d]==0){
                neighbour_x[d] += cell.box[d];
            }
        }

        bool cached=false;

        unsigned cached_bond_index=cached_bonds.size(); // This is the index it will have, _if_ it is cached
        for(unsigned bead_i=0; bead_i < resident.size(); bead_i++){
            auto &bead=resident[bead_i];

            //fprintf(stderr, "    Recv : bead_i=%u, home=%u, other=%u\n", bead_i, get_hash_code(bead.id), get_hash_code(incoming.id));

            // This implicitly interacts each bead with itself, which is handled with a
            // distance check in calc_force.
            float dx[3];
            vec3_sub(dx, bead.x, neighbour_x);
            float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            if(dr_sqr >=1 || dr_sqr < MIN_DISTANCE_CUTOFF_SQR){ // The min threshold avoid large forces, and also skips self-interaction
                continue;
            }
            float dr=pow_half(dr_sqr);

            float kappa=0.0f;
            float r0=cell.bond_r0;
            if(!(is_monomer(incoming.id) || is_monomer(bead.id))){
                if(get_polymer_id(bead.id) == get_polymer_id(incoming.id) ){
                    auto other_polymer_offset=get_polymer_offset(incoming.id);
                    for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
                        if(other_polymer_offset==bead.bond_partners[i]){
                            kappa=cell.bond_kappa;
                            break;
                        }
                    }
                }
            }

            #ifndef NDEBUG
            if(kappa>0){
                cell.debug_hookean_bonds_seen[bead_i]+=1;
            }
            #endif

            auto strength=cell.interactions[MAX_BEAD_TYPES*get_bead_type(bead.id)+get_bead_type(incoming.id)];

            float f[3];
            dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float[3]>(
                cell.inv_root_dt,
                cell.t_hash,
                dx, dr,
                kappa, r0, 
                strength.conservative,
                strength.sqrt_dissipative,
                get_hash_code(bead.id), get_hash_code(incoming.id),
                bead.v, incoming.v,
                f
            );

            vec3_add(bead.f, f);

            if(kappa!=0.0f){
                //std::cerr<<"K: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
                if(!cached){
                    raw_cached_bond_t tmp;
                    tmp.bead_hash=get_hash_code(incoming.id).hash;
                    vec3_copy(tmp.x, neighbour_x);
                    cached_bonds.push_back(tmp);
                    //std::cerr<<" caching, bead_i="<<bead_i<<", bead_id="<<get_hash_code(bead.id)<<" target="<<get_hash_code(incoming.id)<<"\n";
                    cached=true;
                }

                static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
                bool hit = get_polymer_offset(incoming.id) == bead.angle_bonds[0].partner_head || get_polymer_offset(incoming.id) == bead.angle_bonds[0].partner_tail;
                if(hit){
                    //std::cerr<<"  CB: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
                    if(cell.cached_bond_indices[bead_i]==0xFF){
                        //std::cerr<<"  S\n";
                        cell.cached_bond_indices[bead_i] = cached_bond_index;
                    }else{
                        #ifndef NDEBUG
                        assert(cell.cached_bond_indices[bead_i] != 0xFE); // assert angle not already done
                        #endif

                        // Once both partners have arrived, we calculate force and push onto outgoing_forces
                        //std::cerr<<"Force!`\n";
                        assert(bead.id == resident[bead_i].id);
                        calc_angle_force_for_middle_bead<EnableLogging>(cell, bead, cell.cached_bond_indices[bead_i], cached_bond_index);
                        assert(force_outgoing.size()>0);
                        cell.rts |= RTS_FLAG_force;
                        //std::cerr<<"DoenFr\n";
                        #ifndef NDEBUG
                        cell.cached_bond_indices[bead_i] = 0xFE; // Mark angle as done
                        #endif
                    }
                }
            }
        }
    }

    template<bool EnableLogging>
    static void calc_angle_force_for_middle_bead(device_state_t &cell, raw_bead_resident_t &bead, unsigned cache_index_a, unsigned cache_index_b)
    {
        assert(cell.phase==SharingAndForcing);

        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);

        const int ai=0;
        // We only call this function if we already know it is an angle bond
        assert(bead.angle_bonds[ai].partner_head!=0xFF);

        BeadHash head_hash=make_hash_from_offset(bead.id, bead.angle_bonds[ai].partner_head);
        BeadHash tail_hash=make_hash_from_offset(bead.id, bead.angle_bonds[ai].partner_tail);

        const auto *head=&cached_bonds[cache_index_a];
        const auto *tail=&cached_bonds[cache_index_b];

        if( BeadHash{head->bead_hash}.reduced_equals(tail_hash) ){
            std::swap(head, tail);
        }
        
        // The cache copies should already have wrapping applied
        float first[3], second[3];
        vec3_sub(first, bead.x, head->x);
        vec3_sub(second, tail->x, bead.x);

        float FirstLength= vec3_l2_norm(first);
        float SecondLength  = vec3_l2_norm(second);
        assert(FirstLength < 1);
        assert(SecondLength < 1);

        float headForce[3], middleForce[3], tailForce[3];

        dpd_maths_core_half_step_raw::calc_angle_force<EnableLogging,float,float[3],float[3]>(
            (float)bead.angle_bonds[ai].kappa, 0.0f, 0.0f,
            first, FirstLength,
            second, SecondLength,
            headForce, middleForce, tailForce
        );

        vec3_add(bead.f, middleForce);

        force_outgoing.alloc_back();
        force_outgoing.back().target_hash=head_hash.hash;
        vec3_copy(force_outgoing.back().f, headForce);

        force_outgoing.alloc_back();
        force_outgoing.back().target_hash=tail_hash.hash;
        vec3_copy(force_outgoing.back().f, tailForce);

        if(EnableLogging && ForceLogging::logger()){
            ForceLogging::logger()->LogBeadTripleProperty(BeadHash{head->bead_hash}, BeadHash{bead.id}, BeadHash{tail->bead_hash}, "f_next_angle_head", 3, headForce);
            ForceLogging::logger()->LogBeadTripleProperty(BeadHash{head->bead_hash}, BeadHash{bead.id}, BeadHash{tail->bead_hash}, "f_next_angle_mid", 3, middleForce);
            ForceLogging::logger()->LogBeadTripleProperty(BeadHash{head->bead_hash}, BeadHash{bead.id}, BeadHash{tail->bead_hash}, "f_next_angle_tail", 3, tailForce);
        }
    }


    static void on_send_force(device_state_t &cell,raw_force_input_t &outgoing)
    {
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

        assert(!force_outgoing.empty());
        copy_force_input(&outgoing, &force_outgoing.back());
        force_outgoing.pop_back();

        if(force_outgoing.empty()){
            cell.rts &= ~RTS_FLAG_force;
        }
    }

    static void on_recv_force(device_state_t &cell, const raw_force_input_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        for(auto &b : resident){
            if( BeadHash{incoming.target_hash}.reduced_equals(BeadHash{b.id})){
                vec3_add(b.f, incoming.f);
                if(ForceLogging::logger()){
                    ForceLogging::logger()->LogBeadProperty(BeadHash{b.id}, "recv_force_angle", 3, incoming.f);
                }

            }
        }
    }

};

#endif
