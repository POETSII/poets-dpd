#ifndef basic_dpd_engine_v5_raw_handlers_hpp
#define basic_dpd_engine_v5_raw_handlers_hpp

// Add orchestrator specific-patch to allow values in handler_log
// No effect for other platforms.
#ifdef TINSEL
#include "mini_printf.hpp"
#endif

#include "dpd/storage/bag_wrapper.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"
#include "dpd/core/edge_wrap.hpp"

#include <iterator>
#include <cstring>
#include <fstream>
#include <array>

// Work around lack of c++17 in Orchestrator
template <class T, std::size_t N>
constexpr std::size_t std_size(const T (&array)[N]) noexcept
{
    return N;
}

inline void tinsel_require(bool cond, const char *msg)
{
    if(!cond){
        #ifdef PDPD_TINSEL
        puts(msg);
        #else
        throw std::runtime_error(msg);
        #endif
    }
}

template<bool NoBonds=false,bool NoRandom=false>
struct BasicDPDEngineV5RawHandlersImpl
{
    static constexpr const char *THIS_HEADER=__FILE__;

    static constexpr size_t MAX_BONDS_PER_BEAD = 4;
    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_ANGLE_BONDS_PER_BEAD=1;
    static constexpr size_t MAX_CACHED_BONDS_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = MAX_BEADS_PER_CELL * 3; // TODO : This seems very pessimistic

    static constexpr size_t MAX_BEAD_TYPES=8;

    static constexpr bool NO_BONDS=NoBonds;

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winvalid-offsetof"
    // TODO : g++ correctly complains, as BasicDPDEngine::bead_resident_t is not strictly POD

    template<class A>
    static constexpr bool is_bead_view()
    {
        constexpr size_t OidA=offsetof(A,id);
        return
            ((OidA%4)==0)
            && (offsetof(A,x)==OidA+4)
            && (offsetof(A,v)==OidA+16)
            && (sizeof(A::v)==12);
    }


    template<class A>
    static constexpr bool is_bead_resident()
    {
        constexpr size_t OidA=offsetof(A,id);
        return is_bead_view<A>()
            && (offsetof(A,bond_partners)==OidA+10*4)
            && (offsetof(A,angle_bonds)==OidA+10*4+4)
            && (offsetof(A,t)==offsetof(A,angle_bonds)+4)
            && (offsetof(A,checksum)==offsetof(A,t)+4
            && (sizeof(A::checksum)==4)
        );
    }
    
   template<class A>
    static constexpr bool is_force_input()
    {
        constexpr size_t O=offsetof(A,target_hash);
        return ((O%4)==0) && (offsetof(A,f)==O+4) && (sizeof(A::f)==12);
    }

    struct copy_bead_view
    {
        template<class A, class B>
        static void copy(A *dst, const B *src)
        {
            static_assert(is_bead_view<A>());
            static_assert(is_bead_view<B>());

            memcpy32<(1+3+3)>(&dst->id, &src->id);
            //memcpy32(&dst->id, &src->id,(1+3+3));
        }
    };



    struct copy_bead_resident
    {
        template<class A, class B>
        static void copy(A *dst, const B *src)
        {
            static_assert(is_bead_resident<A>());
            static_assert(is_bead_resident<B>());
            static_assert(offsetof(A,checksum)-offsetof(A,id)==13*4);
            memcpy32<(1+3+3+3+4)>((uint32_t*)&dst->id, (uint32_t*)&src->id);
            //memcpy32((uint32_t*)&dst->id, (uint32_t*)&src->id, (1+3+3+3+4));
        }
    };

    #pragma GCC diagnostic pop


    struct copy_cached_bond
    {
        template<class A, class B>
        static void copy(A *dst, const B *src)
        {
            /*dst->id=src->id;
            for(int i=0; i<3; i++){
                dst->x[i]=src->x[i];
            }
            */
           memcpy32<4>((uint32_t*)dst, (const uint32_t*)src);
        }
    };


    struct raw_bead_view_t
    {
        uint32_t id;
        float x[3];
        float v[3];

        template<class T>
        void operator=(const T &x)
        {
            copy_bead_view(this, &x);
        }

        template<class Visitor>
        void walk(Visitor &vv)
        {
            vv("id", id);
            vv("x", x, std_size(x));
            vv("v", v, std_size(v));
        };
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

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("partner_head", partner_head);
                vv("partner_tail", partner_tail);
                vv("kappa", kappa);
                vv("_pad_", _pad_);
            }
        }angle_bonds[MAX_ANGLE_BONDS_PER_BEAD];

        uint32_t t;  // Time-step
        uint32_t checksum;

        void operator=(const raw_bead_resident_t &x)
        {
            copy_bead_resident::copy(this, &x);
        }

        template<class Visitor>
        void walk(Visitor &vv)
        {
            vv("id", id);
            vv("x", x, std_size(x));
            vv("v", v, std_size(v));
            vv("f", f, std_size(f));
            vv("bond_partners", bond_partners, std_size(bond_partners));
            vv("angle_bonds", angle_bonds, std_size(angle_bonds));
            vv("t", t);
            vv("checksum", checksum);
        }
    };
    static_assert(sizeof(raw_bead_resident_t)==56);

    using raw_angle_bond_info_t = decltype(raw_bead_resident_t::angle_bonds[0]);


    struct raw_cached_bond_t
    {
        uint32_t bead_hash;
        float x[3];

        template<class Visitor>
        void walk(Visitor &vv)
        {
            vv("bead_hash", bead_hash);
            vv("x", x, std_size(x));
        }
    };


    struct raw_force_input_t
    {
        uint32_t target_hash;
        float f[3];

        void operator=(const raw_force_input_t &x)
        {
            memcpy32<sizeof(*this)/4>((uint32_t*)this, (const uint32_t*)&x);
            //memcpy32((uint32_t*)this, (const uint32_t*)&x, sizeof(*this)/4);
        }

        template<class Visitor>
        void walk(Visitor &vv)
        {
            vv("target_hash", target_hash);
            vv("f", f, std_size(f));
        }
    };

    template<class A, class B>
    static void copy_force_input(A *dst, const B *src)
    {
        static_assert(is_force_input<A>());
        static_assert(is_force_input<B>());
        memcpy32<sizeof(B)/4>(&dst->target_hash, &src->target_hash);
        //memcpy32(&dst->target_hash, &src->target_hash, sizeof(B)/4);
    }

    struct device_state_t
    {
        float box[3];
        int32_t location[3];
        uint32_t edge_bits;
        float dt;
        float scaled_inv_root_dt;
        float bond_r0;
        float bond_kappa;
        struct {
            float conservative;
            float sqrt_dissipative;

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("conservative", conservative);
                vv("sqrt_dissipative", sqrt_dissipative);
            }
        }interactions[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
        uint32_t t;
        uint64_t t_hash;
        uint64_t t_seed;
        uint32_t interval_size;
        uint32_t output_reps;  // Used to get around lost output packets. Each packet sent this many times
        

        int32_t phase;  // of type Phase
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

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("elements", elements, std_size(elements));
                vv("n", n);
                vv("lost", lost);
            }
        }resident;
        struct{
            raw_bead_resident_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("elements", elements, std_size(elements));
                vv("n", n);
                vv("lost", lost);
            }
        }migrate_outgoing;
        struct force_input_bag
        {
            raw_force_input_t elements[MAX_OUTGOING_FORCES_PER_CELL];
            uint16_t n;
            uint16_t lost;

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("elements", elements, std_size(elements));
                vv("n", n);
                vv("lost", lost);
            }
        } force_outgoing;
        struct cached_bond_bag
        {
            raw_cached_bond_t elements[MAX_CACHED_BONDS_PER_CELL];
            uint16_t n;
            uint16_t lost;

            template<class Visitor>
            void walk(Visitor &vv)
            {
                vv("elements", elements, std_size(elements));
                vv("n", n);
                vv("lost", lost);
            }
        } cached_bonds;

        // Should probably be co-located with the beads for caching
        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
        uint8_t cached_bond_indices[MAX_BEADS_PER_CELL]; // 0xff if not seen, otherwise the first bond that was seen
    
        template<class Visitor>
        void walk(Visitor &vv)
        {
            vv("box", box, std_size(box));
            vv("location", location, std_size(location));
            vv("edge_bits", edge_bits);
            vv("dt", dt);
            vv("scaled_inv_root_dt", scaled_inv_root_dt);
            vv("bond_r0", bond_r0);
            vv("bond_kappa", bond_kappa);
            vv("interactions", interactions, std_size(interactions));
            vv("t", t);
            vv("t_hash", t_hash);
            vv("t_seed", t_seed);
            vv("interval_size", interval_size);
            vv("output_reps", output_reps);
            
            vv("phase", phase);
            vv("rts", rts);
            vv("intervals_todo", intervals_todo);
            vv("interval_offset", interval_offset);
            vv("output_reps_todo", output_reps_todo);
            vv("outputs_todo", outputs_todo);
            vv("share_todo", share_todo);
            vv("resident", resident);
            vv("migrate_outgoing", migrate_outgoing);
            vv("force_outgoing", force_outgoing);
            vv("cached_bonds", cached_bonds);
            vv("cached_bond_indices", cached_bond_indices, std_size(cached_bond_indices));
        }
    };

    static const bool DO_CHECKSUM=false;

    template<class TDeviceState=device_state_t>
    static uint32_t calc_rts(const TDeviceState &cell)
    {
        return cell.rts;
    }

    template<class TDeviceState=device_state_t>
    static bool on_barrier(TDeviceState &cell)
    {
        switch(cell.phase){
            default: assert(false); // fallthrough
            case PreMigrate:
            case Outputting: 
            case SharingAndForcing: return on_barrier_pre_migrate(cell); break;
            case Migrating: return on_barrier_pre_share(cell); break;            
        }
    }

    template<class TDeviceState=device_state_t>
    static void /*__attribute__ ((noinline)) __attribute__((optimize("O0")))*/ on_init(TDeviceState &cell)
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


    template<class TDeviceState=device_state_t>
    static bool on_barrier_pre_migrate(TDeviceState &cell)
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
                migrate_outgoing.template push_back_copier<copy_bead_resident>(b);

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

    template<class TDeviceState=device_state_t, class TRawBeadResident=raw_bead_resident_t>
    static void on_send_migrate(TDeviceState &cell, TRawBeadResident &outgoing)
    {
        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);
        assert(!migrate_outgoing.empty());

        migrate_outgoing.template pop_back<copy_bead_resident>(outgoing);

        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;
    }

    template<class TDeviceState=device_state_t, class TRawBeadResident=raw_bead_resident_t>
    static void on_recv_migrate(TDeviceState &cell, const TRawBeadResident &incoming)
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
            //printf("  recv_migrate : accept.\n");
            resident.template push_back_copier<copy_bead_resident>(incoming);
        }else{
            //printf("  recv_migrate, loc=[%d,%d,%d], point=[%d,%d,%d] : reject.\n", cell.location[0], cell.location[1], cell.location[2], incoming_loc[0], incoming_loc[1], incoming_loc[2]);
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

    template<class TDeviceState=device_state_t>
    static bool on_barrier_pre_share(TDeviceState &cell)
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

    template<class TDeviceState=device_state_t, class TRawBeadView=raw_bead_view_t>
    static void on_send_share(TDeviceState &cell, TRawBeadView &outgoing)
    {
        auto resident=make_bag_wrapper(cell.resident);

        assert(cell.share_todo>0);
        --cell.share_todo;
        const auto &b = resident[cell.share_todo];

        copy_bead_view::copy(&outgoing, &b);
        outgoing.id = b.id;

        if(cell.share_todo==0){
            cell.rts &= ~RTS_FLAG_share;
        }
    }

    template<bool EnableLogging=false, class TDeviceState=device_state_t, class TRawBeadView=raw_bead_view_t>
    static void on_recv_share(TDeviceState &cell, const TRawBeadView &incoming)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<"), p="<<&cell<<", nres="<<cell.resident.n<<", other="<<get_hash_code(incoming.id)<<"\n";

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        
        float neighbour_x[3]={
            incoming.x[0],
            incoming.x[1],
            incoming.x[2]
        };
        do_neighbour_wrap(neighbour_x, cell.edge_bits, cell.box);

        bool cached=false;

        const float ONE=1.0f;
        const float R0=cell.bond_r0; 
        const float KAPPA=cell.bond_kappa;

        unsigned cached_bond_index;
        if constexpr(!NO_BONDS){
            cached_bond_index=cached_bonds.size(); // This is the index it will have, _if_ it is cached
        }
        for(unsigned bead_i=0; bead_i < resident.size(); bead_i++){
            auto &bead=resident[bead_i];

            //fprintf(stderr, "    Recv : bead_i=%u, home=%u, other=%u\n", bead_i, get_hash_code(bead.id), get_hash_code(incoming.id));

            // This implicitly interacts each bead with itself, which is handled with a
            // distance check in calc_force.
            float dx[3];
            float dr_sqr;
            if(0){
                vec3_sub(dx, bead.x, neighbour_x);
                dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr_sqr >=1 || dr_sqr < MIN_DISTANCE_CUTOFF_SQR){ // The min threshold avoid large forces
                    continue;
                }
            }else{
                // Marginally faster
                dx[0] = bead.x[0] - neighbour_x[0];
                dr_sqr = dx[0]*dx[0];
                if(dr_sqr >= ONE){
                    continue;
                }
                dx[1] = bead.x[1] - neighbour_x[1];
                dr_sqr += dx[1]*dx[1];
                if(dr_sqr >= ONE){
                    continue;
                }
                dx[2] = bead.x[2] - neighbour_x[2];
                dr_sqr += dx[2]*dx[2];
                if(dr_sqr >= ONE){
                    continue;
                }
            }
            if(bead.id==incoming.id){
                continue;
            }
            float dr=pow_half(dr_sqr);

            bool is_bonded=false;
            if constexpr(!NO_BONDS){
                if(!(is_monomer(incoming.id) || is_monomer(bead.id))){
                    if(get_polymer_id(bead.id) == get_polymer_id(incoming.id) ){
                        auto other_polymer_offset=get_polymer_offset(incoming.id);
                        for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
                            if(other_polymer_offset==bead.bond_partners[i]){
                                is_bonded=true;
                                break;
                            }
                        }
                    }
                }
            }

            auto bead_type1=get_bead_type(bead.id);
            auto bead_type2=get_bead_type(incoming.id);
            auto strength=cell.interactions[MAX_BEAD_TYPES*bead_type1+bead_type2];

            float f[3];
            const auto MathOptions = NoRandom
                ? dpd_maths_core_half_step_raw::Options::DisableRandom
                : dpd_maths_core_half_step_raw::Options::Default;
            dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float[3],MathOptions>(
                cell.scaled_inv_root_dt,
                cell.t_hash,
                dx, dr,
                is_bonded, KAPPA, R0, 
                strength.conservative,
                strength.sqrt_dissipative,
                BeadHash{bead.id}, BeadHash{incoming.id},
                bead.v, incoming.v,
                f
            );

            vec3_add(bead.f, f);

            if constexpr(!NO_BONDS){
                if(is_bonded){
                    //std::cerr<<"K: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
                    if(!cached){
                        auto &tmp=cached_bonds.alloc_back();
                        tmp.bead_hash=get_hash_code(incoming.id).hash;
                        vec3_copy(tmp.x, neighbour_x);
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
    }

    template<class TDeviceState=device_state_t, class TRawBeadResident=raw_bead_resident_t>
    static void calc_angle_force_for_middle_bead(TDeviceState &cell, TRawBeadResident &bead, unsigned cache_index_a, unsigned cache_index_b)
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


    template<class TDeviceState=device_state_t, class TRawForceInput=raw_force_input_t>
    static void on_send_force(TDeviceState &cell,TRawForceInput &outgoing)
    {
        if constexpr(NO_BONDS){
            assert(0);
        }else{
            auto force_outgoing=make_bag_wrapper(cell.force_outgoing);

            assert(!force_outgoing.empty());
            copy_force_input(&outgoing, &force_outgoing.back());
            force_outgoing.pop_back();

            if(force_outgoing.empty()){
                cell.rts &= ~RTS_FLAG_force;
            }
        }
    }

    template<class TDeviceState=device_state_t, class TRawForceInput=raw_force_input_t>
    static void on_recv_force(TDeviceState &cell, const TRawForceInput &incoming)
    {
        if constexpr(NO_BONDS){
            assert(0);
        }else{
            auto resident=make_bag_wrapper(cell.resident);

            for(auto &b : resident){
                if( BeadHash{b.id}.reduced_equals( BeadHash{incoming.target_hash})){
                    vec3_add(b.f, incoming.f);
                }
            }
        }
    }

    template<class TDeviceState=device_state_t, class TRawBeadResident=raw_bead_resident_t>
    static void on_send_output(TDeviceState &cell, TRawBeadResident &outgoing)
    {
        assert(cell.phase == Outputting);
        assert(cell.interval_offset==cell.interval_size && cell.outputs_todo>0);
        assert(cell.output_reps_todo>0);
        assert(cell.force_outgoing.n==0);
        
        cell.outputs_todo--;
        copy_bead_resident::copy(&outgoing, cell.resident.elements+cell.outputs_todo);

        if(cell.outputs_todo==0){
            cell.output_reps_todo -= 1;
            if(cell.output_reps_todo > 0){
                cell.outputs_todo=cell.resident.n;
            }
        }
        
        cell.rts=cell.outputs_todo==0 ? 0 : RTS_FLAG_output;
    }

};

using BasicDPDEngineV5RawHandlers = BasicDPDEngineV5RawHandlersImpl<false>;
using BasicDPDEngineV5RawHandlersNoBonds = BasicDPDEngineV5RawHandlersImpl<true>;

#endif
