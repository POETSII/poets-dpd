#ifndef hemi_dpd_engine_v1_raw_handlers_hpp
#define hemi_dpd_engine_v1_raw_handlers_hpp

// Add orchestrator specific-patch to allow values in handler_log
// No effect for other platforms.
//#include "mini_printf.hpp"

#include "dpd/storage/bag_wrapper.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"
#include "dpd/core/edge_wrap.hpp"

#include <iterator>
#include <cstring>
#include <fstream>
#include <array>


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

struct HemiDPDEngineV1RawHandlers
{
    static constexpr const char *THIS_HEADER=__FILE__;

    static constexpr size_t MAX_BEADS_PER_CELL = 32;
    static constexpr size_t MAX_BEAD_TYPES=8;
    static constexpr size_t MAX_OUTGOING_FORCES_PER_CELL = 27 * MAX_BEADS_PER_CELL;

    static constexpr size_t MAX_BEADS_PER_MESSAGE = 1;
    static constexpr size_t MAX_VIEWS_PER_MESSAGE = 2;
    static constexpr size_t MAX_FORCES_PER_MESSAGE = 3;

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
    struct bead_view_t
    {
        uint32_t id;
        float x[3];
        float v[3];

        void operator=(const bead_view_t &x)
        {
            memcpy32(*this, x);
        }
    };
    static_assert(sizeof(bead_view_t)==28);

    struct bead_resident_t
        : public bead_view_t
    {
        // Extras
        float f[3];

        uint32_t t;  // Time-step

        void operator=(const bead_resident_t &x)
        {
            memcpy32(*this, x);
        }
    };
    static_assert(sizeof(bead_resident_t)==44);

    struct force_input_t
    {
        uint32_t target_hash;
        float f[3];

        void operator=(const force_input_t &x)
        {
            memcpy32(*this, x);
        }
    };
    static_assert(sizeof(force_input_t)==16);

    struct message_t
    {
        uint16_t type;
        uint16_t n;
        union {
            bead_resident_t beads[MAX_BEADS_PER_MESSAGE];
            bead_view_t views[MAX_VIEWS_PER_MESSAGE];
            force_input_t forces[MAX_FORCES_PER_MESSAGE];
        };
    };
    static_assert(sizeof(message_t) <= 60);

    struct device_state_t
    {
        float box[3];
        int32_t location[3];
        uint32_t edge_bits;
        float dt;
        float scaled_inv_root_dt;
        struct {
            float conservative;
            float sqrt_dissipative;
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
            bead_resident_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }resident;
        struct{
            bead_resident_t elements[MAX_BEADS_PER_CELL];
            uint16_t n;
            uint16_t lost;
        }migrate_outgoing;
        struct force_input_bag
        {
            force_input_t elements[MAX_OUTGOING_FORCES_PER_CELL];
            uint16_t n;
            uint16_t lost;
        } force_outgoing;
    };

    template<class device_state_t=device_state_t>
    static uint32_t calc_rts(const device_state_t &cell)
    {
        if(cell.rts==0){
            assert(cell.force_outgoing.n==0);
        }
        return cell.rts;
    }

    template<class device_state_t=device_state_t>
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

    template<class device_state_t=device_state_t>
    static void /*__attribute__ ((noinline)) __attribute__((optimize("O0")))*/ on_init(device_state_t &cell)
    {
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);
    }

    static bool on_barrier_pre_migrate(device_state_t &cell)
    {
        //fprintf(stderr, "b_pre_migrate, t=%u, loc=[%d,%d,%d]\n", cell.t, cell.location[0], cell.location[1],cell.location[2]);

        assert(cell.phase==PreMigrate || cell.phase==SharingAndForcing || cell.phase==Outputting);

        //int tmp;
        //printf("on_barrier_pre_migrate : %x, sp=%x\n", (unsigned)(intptr_t)&cell, (unsigned)(intptr_t)&tmp);

        assert(cell.rts==0);
        assert(make_bag_wrapper(cell.force_outgoing).empty());

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

        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        assert(migrate_outgoing.empty());

        int i=resident.size();
        while(0 < i){
            --i;
            bead_resident_t &b = resident[i];

            // We always apply the mom step here, rather than doing it in 
            // a seperate step at the end. This means that the input needs
            // to have been "reverse" update_mom'd, and the output must be
            // mom'd again.
            dpd_maths_core_half_step_raw::update_mom(cell.dt, b);

            // Actual position update, clears force
            dpd_maths_core_half_step_raw::update_pos(cell.dt, cell.box, b);

            ++b.t;
        
            // Check it is still in the right cell
            int32_t true_loc[3];
            vec3_floor_nn(true_loc, b.x);
            if( !vec3_equal(true_loc, cell.location) ){
                //fprintf(stderr, "Bead %u leaving [%d,%d,%d] for [%d,%d,%d]\n", b.id, cell.location[0], cell.location[1], cell.location[2],  true_loc[0], true_loc[1], true_loc[2]);

                // Movement is fairly unlikely
                memcpy32( migrate_outgoing.alloc_back(), b );

                memswap32(resident.back(), b); // A waste if this is the last bead, but makes code simpler and smaller
                resident.pop_back();
            }
        }

        cell.t_hash = get_t_hash(cell.t, cell.t_seed);
        cell.t += 1;

        cell.phase=Migrating;
        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;

        //fprintf(stderr, "   migrate.n=%d, rts=%x\n", cell.migrate_outgoing.n, cell.rts);

        return true;
    }

    static void on_send_migrate(device_state_t &cell, message_t &outgoing)
    {
        static_assert(MAX_BEADS_PER_MESSAGE==1);

        auto migrate_outgoing=make_bag_wrapper(cell.migrate_outgoing);
        assert(!migrate_outgoing.empty());

        outgoing.type=RTS_INDEX_migrate;
        outgoing.n=1;
        memcpy32(outgoing.beads[0], migrate_outgoing.back());
        migrate_outgoing.pop_back();

        //fprintf(stderr, "Bead %u migrating out of [%d,%d,%d] \n", outgoing.beads[0].id, cell.location[0], cell.location[1], cell.location[2]);

        cell.rts=migrate_outgoing.empty() ? 0 : RTS_FLAG_migrate;
    }

    static void on_recv_migrate(device_state_t &cell, const message_t &incoming)
    {
        static_assert(MAX_BEADS_PER_MESSAGE==1);
        auto &b=incoming.beads[0];

        if(cell.location[0] == floor_nn(b.x[0])){
            if(cell.location[1] == floor_nn(b.x[1])){
                if(cell.location[2] == floor_nn(b.x[2])){
                    //fprintf(stderr, "Bead %u migrating into [%d,%d,%d] at %u \n", b.id, cell.location[0], cell.location[1], cell.location[2], cell.t);

                    
                    for(unsigned i=0; i<cell.resident.n; i++){
                        assert( cell.resident.elements[i].id != b.id );
                    }

                    assert( cell.resident.n < MAX_BEADS_PER_CELL );
                    memcpy32( cell.resident.elements[cell.resident.n], b );
                    cell.resident.n += 1;
                    //assert(cell.resident.n<2);
                }        
            }               
        }
    }

    static bool on_barrier_pre_share(device_state_t &cell)
    {
        //fprintf(stderr, "b_pre_share, t=%u, loc=[%d,%d,%d], migrate.n=%d\n", cell.t, cell.location[0], cell.location[1],cell.location[2], cell.migrate_outgoing.n);

        assert(cell.phase==Migrating);

        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        unsigned nResident=cell.resident.n;

        for(unsigned i=1; i<nResident; i++){
            cell.resident.n=i; // Temporarily pretend there are i residents...
            bead_resident_t &b = cell.resident.elements[i];  // ... and interact/accumulate into the i'th bead
            float *fAcc=b.f;   // We can accumulate directly into the force vector in memory
            interact_bead(cell, b.x[0], b.x[1], b.x[2], b, fAcc);
        }
        cell.resident.n=nResident; // Restore it back
        cell.share_todo = nResident;

        cell.phase=SharingAndForcing;
        cell.rts=cell.share_todo==0 ? 0 : RTS_FLAG_share;

        return true;
    }

    static void on_send_share(device_state_t &cell, message_t &outgoing)
    {
        auto resident=make_bag_wrapper(cell.resident);

        unsigned todo=std::min<unsigned>(cell.share_todo, MAX_VIEWS_PER_MESSAGE);

        const bead_resident_t *src=cell.resident.elements+cell.share_todo-todo;
        for(unsigned i=0; i<todo; i++){
            memcpy32<bead_view_t>( outgoing.views[i], src[i]);
        }


        cell.share_todo -= todo;

        outgoing.type=RTS_INDEX_share;
        outgoing.n = todo;

        if(cell.share_todo==0){
            assert(cell.rts==RTS_FLAG_share);
            if(cell.force_outgoing.n==0){
                cell.rts=0;
            }else{
                cell.rts = RTS_FLAG_force;
            }
        }
    }

    static bool interact_bead(device_state_t &cell, float bx0, float bx1, float bx2, const bead_view_t &incoming, float *fOpposite)
    {

        const float ONE=1.0f;

        bool any_hit=false;

        // Ideally these are in registers
        float fInv[3];

        float dx[3];
        float dr_sqr;

        for(unsigned bead_i=0; bead_i < cell.resident.n; bead_i++){
            auto &bead=cell.resident.elements[bead_i];
            
            // Marginally faster
            dx[0] = bead.x[0] - bx0;
            dr_sqr = dx[0]*dx[0];
            if(dr_sqr >= ONE){
                continue;
            }
            dx[1] = bead.x[1] - bx1;
            dr_sqr += dx[1]*dx[1];
            if(dr_sqr >= ONE){
                continue;
            }
            dx[2] = bead.x[2] - bx2;
            dr_sqr += dx[2]*dx[2];
            if(dr_sqr >= ONE){
                continue;
            }

            float dr=pow_half(dr_sqr);
            if(dr < 0.00001f){
                continue;
            }

            auto bead_type1=BeadHash{bead.id}.get_bead_type();
            auto bead_type2=BeadHash{incoming.id}.get_bead_type();
            auto strength=cell.interactions[MAX_BEAD_TYPES*bead_type1+bead_type2];

            float f[3];
            const bool EnableLogging=false;
            dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float[3]>(
                cell.scaled_inv_root_dt,
                cell.t_hash,
                dx, dr,
                false, 0.0f, 0.0f, 
                strength.conservative,
                strength.sqrt_dissipative,
                BeadHash{bead.id}, BeadHash{incoming.id},
                bead.v, incoming.v,
                f
            );

            vec3_add(bead.f, f);
            if(any_hit){
                vec3_sub(fInv, f);
            }else{
                any_hit=true;
                vec3_neg(fInv, f);
            }
        }

        if(any_hit){
            vec3_copy(fOpposite, fInv );
        }
        return any_hit;
    }

    static void on_recv_share(device_state_t &cell, const message_t &incoming)
    {
        bool hit=false;
        unsigned incoming_n = MAX_VIEWS_PER_MESSAGE==1 ? 1 : incoming.n; 
        for(unsigned i=0; i<incoming_n; i++){        

            // Ideally get promoted into registers
            float nx[3]={
                incoming.views[i].x[0],
                incoming.views[i].x[1],
                incoming.views[i].x[2]
            };

            if(cell.edge_bits){
                do_neighbour_wrap(nx, cell.edge_bits, cell.box);
            }

            auto &f=cell.force_outgoing.elements[cell.force_outgoing.n];
            // Hopefully cell and x vals are passed by reg
            if(interact_bead(cell, nx[0], nx[1], nx[2], incoming.views[i], f.f)){
                f.target_hash = incoming.views[i].id;
                cell.force_outgoing.n += 1;
                hit=true;
            }
        }
        if(hit){
            if(cell.rts==0){ // prioritise share over force
                cell.rts = RTS_FLAG_force;
            }
        }
    }


    static void on_send_force(device_state_t &cell, message_t &outgoing)
    {
        assert(cell.rts==RTS_FLAG_force);
        assert(cell.force_outgoing.n!=0);
        unsigned todo=std::min<unsigned>( cell.force_outgoing.n, std::size(outgoing.forces) );

        auto *dst=&outgoing.forces[0];
        const auto *src=&cell.force_outgoing.elements[cell.force_outgoing.n-todo];

        for(unsigned i=0; i<todo; i++){
            memcpy32<force_input_t>(dst[i], src[i]);
        }
        for(unsigned i=todo; i<std::size(outgoing.forces); i++){
            outgoing.forces[i].target_hash=0xFFFFFFFF;
        }

        cell.force_outgoing.n -= todo;

        outgoing.type=RTS_INDEX_force;
        outgoing.n=todo;

        if(cell.force_outgoing.n==0){
            assert(cell.rts==RTS_FLAG_force);
            cell.rts = 0;
        }
    }

    static void on_recv_force(device_state_t &cell, const message_t &incoming)
    {
        // Try to promote into registers
        uint32_t rs[MAX_FORCES_PER_MESSAGE];
        #pragma GCC unroll 3
        for(unsigned i=0; i<MAX_FORCES_PER_MESSAGE; i++){
            rs[i]=incoming.forces[i].target_hash;
        }

        for(unsigned i=0; i<cell.resident.n; i++){
            uint32_t rid=cell.resident.elements[i].id;

            #pragma GCC unroll 3
            for(unsigned j=0; j<MAX_FORCES_PER_MESSAGE; j++){
                if(rs[j]==rid){
                    vec3_add(cell.resident.elements[i].f, incoming.forces[j].f);
                    break;
                }
            }
        }
    }

    static void on_send_output(device_state_t &cell, message_t &outgoing)
    {
        assert(cell.rts & RTS_FLAG_output);
        assert(cell.phase == Outputting);
        assert(cell.interval_offset==cell.interval_size && cell.outputs_todo>0);
        assert(cell.output_reps_todo>0);
        assert(cell.force_outgoing.n==0);
        
        cell.outputs_todo--;
        memcpy32( outgoing.beads[0], cell.resident.elements[cell.outputs_todo] );
        outgoing.type=RTS_INDEX_output;
        outgoing.n=1;
        
        cell.rts=cell.outputs_todo==0 ? 0 : RTS_FLAG_output;
    }

};


#endif
