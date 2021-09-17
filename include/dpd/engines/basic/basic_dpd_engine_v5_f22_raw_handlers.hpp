#ifndef basic_dpd_engine_v5_f22_raw_handlers_hpp
#define basic_dpd_engine_v5_f22_raw_handlers_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw_handlers.hpp"

struct BasicDPDEngineV5F22RawHandlers
    : BasicDPDEngineV5RawHandlers
{   

    
    
    template<class A, class B>
    static void copy_bead_view(A *dst, const B *src)
    {
        static_assert(offsetof(A,id)==0);
        static_assert(offsetof(A,id)==offsetof(B,id));
        static_assert(offsetof(A,x_f22)==sizeof(A::id));
        static_assert(offsetof(A,x_f22)==offsetof(B,x_f22));
        static_assert(offsetof(A,v)==offsetof(A,x_f22)+sizeof(A::x_f22));
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
        static_assert(offsetof(A,x_f22)==offsetof(B,x_f22));
        static_assert(offsetof(A,v)==offsetof(B,v));
        static_assert(offsetof(A,f)==offsetof(B,f));
        static_assert(offsetof(A,bond_partners)==offsetof(B,bond_partners));
        static_assert(offsetof(A,angle_bonds)==offsetof(B,angle_bonds));
        static_assert(offsetof(A,t)==offsetof(B,t));
        static_assert(offsetof(A,checksum)==offsetof(B,checksum));
        static_assert(sizeof(A)%4==0);
        memcpy32((uint32_t*)dst, (const uint32_t*)src, sizeof(A)/4);
    }

    template<class A, class B>
    static void copy_bead_resident_plain_to_f22(A *dst, const B *src)
    {
        static_assert(sizeof(A)==sizeof(B));
        static_assert(offsetof(A,id)==offsetof(B,id));
        static_assert(offsetof(A,x_f22)==offsetof(B,x));
        static_assert(offsetof(A,v)==offsetof(B,v));
        static_assert(offsetof(A,f)==offsetof(B,f));
        static_assert(offsetof(A,bond_partners)==offsetof(B,bond_partners));
        static_assert(offsetof(A,angle_bonds)==offsetof(B,angle_bonds));
        static_assert(offsetof(A,t)==offsetof(B,t));
        static_assert(offsetof(A,checksum)==offsetof(B,checksum));
        static_assert(sizeof(A)%4==0);
        memcpy32((uint32_t*)dst, (const uint32_t*)src, sizeof(A)/4);
        for(int i=0; i<3; i++){
            dst->x_f22[i] = float_to_f22(src->x[i]);
            //fprintf(stderr, "  dim=%d, in=%f, out=%d\n", i, src->x[i], dst->x_f22[i]);
        }
    }
#pragma GCC diagnostic pop

    struct raw_bead_view_f22_t
    {
        uint32_t id;
        int32_t x_f22[3];
        float v[3];

        template<class T>
        void operator=(const T &x)
        {
            copy_bead_view(this, &x);
        }
    };
    static_assert(sizeof(raw_bead_view_f22_t)==28);

    template<class T>
    static BeadHash get_hash_code(const T &bead)
    { return BeadHash{bead.id}; }

    template<class T>
    static uint32_t get_bead_type(const T &bead)
    { return BasicDPDEngineV5RawHandlers::get_bead_type(bead); }

    

    struct raw_bead_resident_f22_t
    {
        // raw_bead_view_t
        uint32_t id;
        int32_t x_f22[3];
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

        void operator=(const raw_bead_resident_f22_t &x)
        {
            copy_bead_resident(this, &x);
        }
    };
    static_assert(sizeof(raw_bead_resident_f22_t)==56);

    using raw_angle_bond_info_t = decltype(raw_bead_resident_t::angle_bonds[0]);


    struct raw_cached_bond_f22_t
    {
        uint32_t bead_hash;
        int32_t x_f22[3];
    };

    struct device_state_f22_t
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

        int32_t location_f0[3]; // Integer co-ordinates

        Phase phase;
        uint32_t rts;

        uint32_t intervals_todo;
        uint32_t interval_offset; // When offset reaches zero we output

        uint32_t output_reps_todo;
        uint32_t outputs_todo;
        uint32_t share_todo;
        struct{
            raw_bead_resident_f22_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }resident;
        struct{
            raw_bead_resident_f22_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }migrate_outgoing;
        struct force_input_bag
        {
            raw_force_input_t elements[MAX_OUTGOING_FORCES_PER_CELL];
            uint16_t n;
            uint16_t lost;
        } force_outgoing;
        struct cached_bond_bag
        {
            raw_cached_bond_f22_t elements[MAX_CACHED_BONDS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        } cached_bonds;

        // Should probably be co-located with the beads for caching
        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
        uint8_t cached_bond_indices[MAX_BEADS_PER_CELL]; // 0xff if not seen, otherwise the first bond that was seen
    };

    static bool on_barrier(device_state_f22_t &cell)
    {
        switch(cell.phase){
            default: assert(false); // fallthrough
            case PreMigrate:
            case Outputting: 
            case SharingAndForcing: return on_barrier_pre_migrate(cell); break;
            case Migrating: return on_barrier_pre_share(cell); break;            
        }
    }


    static bool on_barrier_pre_migrate(device_state_f22_t &cell)
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

        int32_t box_f22[3];
        for(int d=0; d<3; d++){
            box_f22[d]=cell.box[d]<<22;
        }

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

            dpd_maths_core_half_step_raw::update_pos_f22(cell.dt, box_f22, b);
            //fprintf(stderr, "  nx={%d,%d,%d}\n", b.x_f22[0], b.x_f22[1], b.x_f22[2]);

            ++b.t;
        
            // Check it is still in the right cell
            bool wrong_location=false;
            for(int d=0; d<3; d++){
                wrong_location |= (b.x_f22[d]>>22) != cell.location_f0[d];
            }
            if( wrong_location ){
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

    static void on_recv_migrate(device_state_f22_t &cell, const raw_bead_resident_f22_t &incoming)
    {
        auto resident=make_bag_wrapper(cell.resident);

        bool is_here=true;
        for(int i=0; i<3; i++){
            is_here &= cell.location_f0[i] == (incoming.x_f22[i]>>22);
        }
        if(is_here){
            resident.push_back(incoming);
        }
    }

    static float f22_to_float(int32_t x)
    {
        return float(x) * (1.0/float(1<<22));
    }

    static int32_t float_to_f22(float x)
    {
        return float(x) * float(1<<22);
    }

    template<bool EnableLogging>
    static void on_recv_share(device_state_f22_t &cell, const raw_bead_view_f22_t &incoming)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<"), p="<<&cell<<", nres="<<cell.resident.n<<", other="<<get_hash_code(incoming.id)<<"\n";

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        
        int32_t neighbour_x[3];
        int32_t neighbour_cell_pos[3];
        for(int d=0; d<3; d++){
            neighbour_x[d] = incoming.x_f22[d];

            int neighbour_cell_pos_f0=(neighbour_x[d]>>22);

            if(cell.location_f0[d]==0 && neighbour_cell_pos_f0==cell.box[d]-1){
                neighbour_x[d] -= (cell.box[d]<<22);
            }else if(cell.location_f0[d]==cell.box[d]-1 && neighbour_cell_pos_f0==0){
                neighbour_x[d] += (cell.box[d]<<22);
            }

            neighbour_x[d] += 1<<(22-14-1); // This is a pre-bias for truncation in sqr_f22_to_f28
        }

        bool cached=false;

        unsigned cached_bond_index=cached_bonds.size(); // This is the index it will have, _if_ it is cached
        for(unsigned bead_i=0; bead_i < resident.size(); bead_i++){
            auto &bead=resident[bead_i];

            //fprintf(stderr, "    Recv : bead_i=%u, home=%u, other=%u\n", bead_i, get_hash_code(bead.id), get_hash_code(incoming.id));

            auto sqr_f22_to_f28=[](int32_t x)
            {
                // This is not quite symmetric...
                int32_t tmp=(x >>(22-14));
                return tmp*tmp; // Square to f28
            };

            // This implicitly interacts each bead with itself, which is handled with a
            // distance check in calc_force.
            int32_t dx_f22[3];
            vec3_sub(dx_f22, bead.x_f22, neighbour_x);
            
            int32_t dr_sqr_f28=sqr_f22_to_f28(dx_f22[0]) + sqr_f22_to_f28(dx_f22[1]) + sqr_f22_to_f28(dx_f22[2]);
            if(dr_sqr_f28 >= (1<<28)){
                continue;
            }

            const int32_t MIN_DIST_SQR=1; // We cant use standard threshold, as it rounds to zero
            if(dr_sqr_f28 < MIN_DIST_SQR){ // The min threshold avoid large forces
                continue;
            }
            if(bead.id==incoming.id){
                continue;
            }

            float dr=pow_half( dr_sqr_f28 * float(1.0f/(1<<28)) );
            assert(dr>0);
            float dx[3];
            for(int d=0; d<3; d++){
                dx[d] = f22_to_float(dx_f22[d]);
            }

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
                    raw_cached_bond_f22_t tmp;
                    tmp.bead_hash=get_hash_code(incoming).hash;
                    vec3_copy(tmp.x_f22, neighbour_x);
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

    static void calc_angle_force_for_middle_bead(device_state_f22_t &cell, raw_bead_resident_f22_t &bead, unsigned cache_index_a, unsigned cache_index_b)
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
        for(int d=0; d<3; d++){
            first[d]=f22_to_float(bead.x_f22[d] - head->x_f22[d]);
            second[d]=f22_to_float(tail->x_f22[d] - bead.x_f22[d]);
        }

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


};

#endif
