#ifndef basic_dpd_engine_v5_raw_handlers_hpp
#define basic_dpd_engine_v5_raw_handlers_hpp

#include "dpd/storage/bag_wrapper.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"

#include <iterator>
#include <cstring>

inline void tinsel_require(bool cond, const char *msg)
{
    if(!cond){
        #ifdef TINSEL
        printf(msg);
        #else
        throw std::runtime_error(msg);
        #endif
    }
}

struct BasicDPDEngineV5RawHandlers
{
    static constexpr size_t MAX_BONDS_PER_BEAD = 4;
    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=12;

    static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);

    enum OutputFlags
    {
        RTS_INDEX_migrate=0,
        RTS_INDEX_share=1,
        RTS_INDEX_force=2,
        RTS_INDEX_output=3,

        RTS_FLAG_migrate = 1<<RTS_INDEX_migrate,
        RTS_FLAG_share = 1<<RTS_INDEX_share,
        RTS_FLAG_force = 1<<RTS_INDEX_force,
        RTS_FLAG_output = 1<<RTS_INDEX_output
    };

    enum Phase : uint32_t
    {
        PreMigrate,
        Migrating,
        SharingAndForcing,
        Outputting
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
        return BeadHash{bead_id}.make_reduced_hash_from_polymer_offset(offset);
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

    static BeadHash get_hash_code(const raw_bead_view_t &bead)
    { return BeadHash{bead.id}; }

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
            //This is signed because tinsel can only do int->float conversion. It cannot do unsigned->float
            int8_t kappa;   // Kappa is just an integer. Typically this is quite small, e.g. 5 or 15. Should be less than 127
            uint8_t _pad_;
        }angle_bonds[MAX_ANGLE_BONDS_PER_BEAD];

        uint32_t t;  // Time-step
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
        static_assert( (offsetof(A,v)+sizeof(A::v)) % 4 == 0 );
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
        static_assert(offsetof(A,t)==offsetof(B,t));
        static_assert(offsetof(A,checksum)==offsetof(B,checksum));
        static_assert(sizeof(A)%4==0);
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

    template<class A, class B>
    static void copy_force_input(A *dst, const B *src)
    {
        static_assert(sizeof(A)==sizeof(B));
        static_assert(offsetof(A,target_hash)==offsetof(B,target_hash));
        static_assert( std::is_same<decltype(A::f),decltype(B::f)>::value );
        static_assert( (sizeof(A)%4) == 0 );
        memcpy32((uint32_t*)dst, (const uint32_t*)src, sizeof(B)/4);
    }

    struct device_state_t
    {

        int32_t box[3];
        float dt;
        float inv_root_dt;
        float bond_r0;
        float bond_kappa;
        struct {
            float conservative;
            float sqrt_dissipative;
        }interactions[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
        uint32_t t;
        uint64_t t_hash;
        uint64_t t_seed;
        uint32_t interval_size;
        uint32_t output_reps;  // Used to get around lost output packets. Each packet sent this many times

        int32_t location[3];

        Phase phase;
        uint32_t rts;

        uint32_t intervals_todo;
        uint32_t interval_offset; // When offset reaches zero we output

        uint32_t output_reps_todo;
        uint32_t outputs_todo;
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
        uint8_t cached_bond_indices[MAX_BEADS_PER_CELL]; // 0xff if not seen, otherwise the first bond that was seen
    };

    static const bool DO_CHECKSUM=false;

    static uint32_t calc_rts(const device_state_t &cell)
    {
        return cell.rts;
    }

    static bool on_barrier(device_state_t &cell)
    {
        switch(cell.phase){
            default: assert(false); // fallthrough
            case PreMigrate:
            case Outputting: 
            case SharingAndForcing: return on_barrier_pre_migrate(cell); break;
            case Migrating: return on_barrier_pre_share(cell); break;            
        }
    }

    static void /*__attribute__ ((noinline)) __attribute__((optimize("O0")))*/ on_init(device_state_t &cell)
    {
        if(DO_CHECKSUM){
            //auto resident=make_bag_wrapper(cell.resident);

            uint32_t i=0;
            uint32_t n=cell.resident.n;
            //for(const auto & b : resident){
            while(i<n){
                const auto &b=cell.resident.elements[i];
                if(b.checksum!=calc_checksum(b)){
                    printf("CheckSumInit : nn=%x, n=%x, id=%x, %x, %x\n", (unsigned)n, (unsigned)i, (unsigned)b.id, (unsigned)b.checksum, (unsigned)calc_checksum(b));
                }
                ++i;
            }
        }
    }


    static bool on_barrier_pre_migrate(device_state_t &cell)
    {
        assert(cell.phase==PreMigrate || cell.phase==SharingAndForcing || cell.phase==Outputting);

        //int tmp;
        //printf("on_barrier_pre_migrate : %x, sp=%x\n", (unsigned)(intptr_t)&cell, (unsigned)(intptr_t)&tmp);

        auto resident=make_bag_wrapper(cell.resident);
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);

        if(cell.intervals_todo==0){
            return false;
        }else if(cell.interval_offset==0){
            cell.interval_offset=cell.interval_size;
            --cell.intervals_todo;

            cell.outputs_todo=resident.size();
            cell.phase=Outputting;
            cell.output_reps_todo=cell.output_reps;
            cell.rts=cell.outputs_todo>0 ? RTS_FLAG_output : 0;
            return true;
        }
        cell.interval_offset -= 1;

        make_bag_wrapper(cell.cached_bonds).clear();

        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

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

            ++b.t;

            if(DO_CHECKSUM){
                b.checksum = calc_checksum(b);
            }
        
            // Check it is still in the right cell
            int32_t true_loc[3];
            vec3_floor_nn(true_loc, b.x);
            if( !vec3_equal(true_loc, cell.location) ){
                if(DO_CHECKSUM){
                    for(int i=0; i<3; i++){
                        tinsel_require(0<=true_loc[i] && true_loc[i]<cell.box[i],"BeadOutOfBox");
                        int diff=abs(true_loc[i]-cell.location[i]);
                        tinsel_require( (-1 <= diff && diff<=+1) || (diff==cell.box[i]-1), "BeadJumpedTwoCells");
                    }
                }

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

        return true;
    }


    static void on_send_migrate(device_state_t &cell, raw_bead_resident_t &outgoing)
    {
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);
        assert(!migrate_outgoing.empty());
        
        outgoing=migrate_outgoing.back();
        migrate_outgoing.pop_back();

        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;
    }

    static void on_recv_migrate(device_state_t &cell, const raw_bead_resident_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        if(DO_CHECKSUM){
            if(incoming.checksum != calc_checksum(incoming)){
                printf("CheckSumMigrateRecv\n");
            }
        }

        int32_t incoming_loc[3];
        vec3_floor_nn(incoming_loc, incoming.x);
        if(vec3_equal(incoming_loc , cell.location)){
            resident.push_back(incoming);
        }else{
            if(DO_CHECKSUM){
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
    }

    static bool on_barrier_pre_share(device_state_t &cell)
    {
        assert(cell.phase==Migrating);

        auto resident=make_bag_wrapper(cell.resident);
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        cell.share_todo = resident.size();

        for(unsigned i=0; i<resident.size(); i++){ 
            cell.cached_bond_indices[i]=0xff;
        }

        cell.phase=SharingAndForcing;
        cell.rts=cell.share_todo==0 ? 0 : RTS_FLAG_share;

        return true;
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
            if(dr_sqr >=1 || dr_sqr < MIN_DISTANCE_CUTOFF_SQR){ // The min threshold avoid large forces
                continue;
            }
            if(bead.id==incoming.id){
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

            auto bead_type1=get_bead_type(bead.id);
            auto bead_type2=get_bead_type(incoming.id);
            auto strength=cell.interactions[MAX_BEAD_TYPES*bead_type1+bead_type2];

            float f[3];
            dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float[3]>(
                cell.inv_root_dt,
                cell.t_hash,
                dx, dr,
                kappa, r0, 
                strength.conservative,
                strength.sqrt_dissipative,
                BeadHash{bead.id}, BeadHash{incoming.id},
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
                        // Once both partners have arrived, we calculate force and push onto outgoing_forces
                        //std::cerr<<"Force!`\n";
                        assert(bead.id == resident[bead_i].id);
                        calc_angle_force_for_middle_bead(cell, bead, cell.cached_bond_indices[bead_i], cached_bond_index);
                        cell.rts |= RTS_FLAG_force;
                        //std::cerr<<"DoenFr\n";
                    }
                }
            }
        }
    }

    static void calc_angle_force_for_middle_bead(device_state_t &cell, raw_bead_resident_t &bead, unsigned cache_index_a, unsigned cache_index_b)
    {
        assert(cell.phase==SharingAndForcing);

        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);

        const int ai=0;
        // We only call this function if we already know it is an angle bond
        assert(bead.angle_bonds[ai].partner_head!=0xFF);

        auto head_hash=make_hash_from_offset(bead.id, bead.angle_bonds[ai].partner_head);
        auto tail_hash=make_hash_from_offset(bead.id, bead.angle_bonds[ai].partner_tail);

        const auto *head=&cached_bonds[cache_index_a];
        const auto *tail=&cached_bonds[cache_index_b];

        if( BeadHash{head->bead_hash}.reduced_equals(tail_hash)){
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

        dpd_maths_core_half_step_raw::calc_angle_force<true,float,float[3],float[3]>(
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
            if( BeadHash{b.id}.reduced_equals( BeadHash{incoming.target_hash})){
                vec3_add(b.f, incoming.f);
            }
        }
    }

    static void on_send_output(device_state_t &cell, raw_bead_resident_t &outgoing)
    {
        assert(cell.phase == Outputting);
        assert(cell.interval_offset==cell.interval_size && cell.outputs_todo>0);
        assert(cell.output_reps_todo>0);
        
        cell.outputs_todo--;
        copy_bead_resident(&outgoing, &cell.resident.elements[cell.outputs_todo]);

        if(cell.outputs_todo==0){
            cell.output_reps_todo -= 1;
            if(cell.output_reps_todo > 0){
                cell.outputs_todo=cell.resident.n;
            }
        }
        
        cell.rts=cell.outputs_todo==0 ? 0 : RTS_FLAG_output;
    }

};

#endif
