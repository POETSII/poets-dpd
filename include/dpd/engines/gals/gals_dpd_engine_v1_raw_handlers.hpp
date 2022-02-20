#ifndef gals_dpd_engine_v1_handlers_hpp
#define gals_dpd_engine_v1_handlers_hpp

#include "dpd/storage/bag_wrapper.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"

#include <iterator>
#include <cstring>

#ifndef assert_implies
#define assert_implies(a,b) assert( (a) ? (b) : true )
#endif

/*
    Each device moves through a number of discrete time-steps t=0,1,2,3,...
    Each time-step has an associated round, which contains the working state
    for that round.

    Each device has two active rounds (time-steps):
    - curr = rounds[t] : the time-step the device is currently in
    - next = rounds[t+1] : the time-step the device will move to next
    In principle all the earlier and later rounds "exist", but we only need
    to track these two active rounds (i.e. hold them in memory).

    Each round contains certain sets and properties, which we will define in terms of full beads where convenient, even if only parts are needed:
    - resident : the set of beads which are located within the device for this round
    - migrating : the set of beads which were located in this device for the previous round, but have left for this round
    - shared : beads which have been shared with neighbouring devices
    - viewed : the set of beads received from other devices in this round, not including this device
    - outputs : set of beads at time round.t+1 we need to output before moving to the next round
    - expected_views : a bag of integers describing the expected number of shared beads from neighbours (not including this device)
    - expected_migrations : a bag of integers describing the expected number of migrations from neighbours
    - migrations_received_count : integer giving the number of migration _messages_ which have arrived (not the number that stayed here)
    - have_sent_expected_views : boolean indicating if we have sent the expected views
    - have_sent_expected_migrations : boolean indicating if we have send the expected migrations

    These sets are not fixed within a round:
    - The resident set may grow when beads are migrated into this set
    - The migrating set may shrink as migrating beads are sent to other devices
    - The shared set will grow as resident beads are migrated
    - The viewed set will grow as shared beads are received
    - The outputs set my grow when the round becomes complete, then shrink as outputs are sent
    - expected_views will grow as we receive information from neighbours
    - expected_migrations will grow as we receive information from neighbours 

    In practice information about expected views and migrations will be carried in other
    share and migration messages, but here we considered them as distinct.

    Invariants on all rounds at all times:
    - shared \subsetof resident : We can only share beads that we definitely own in this round
    - resident \intersection migrating = 0 : Beads can only be resident or migrating
        -Note that the migrating beads will eventually come back in the form of shares, but that is handled normally
    - viewed \intersection resident = 0 : Resident beads do not appear in the viewed set
    - outputs \subsetof resident : We only output beads that we own within this round

    A round is considered resident_complete if the following is true:
            |r.expected_migrations| == 26 and    // We must have received full set of expected counts
            sum(r.expected_migrations) == r.migrations_received_count   // every expected migration must have arrived
    
    Once a round is resident complete, the resident set is completely stable and will never change.

    Similarly, we have view_complete(r):
        |r.expected_views| == 26 and
        sum(r.expected_views) == |r.views|

    A round is considered send_complete if the following is true:
        send_complete(r) =
            resident_complete(r) and
            // Our neighbours know how many messages we will send
            r.have_sent_expected_views && r.have_sent_expected_migrations and
            // Make sure the resident set is shared
            r.shared==r.resident and
            r.migrating==0

    A round is considered recv_complete if the following is true:
        recv_complete(r) = 
            resident_complete(r) and
            // Make sure the view set is complete
            |r.expected_views| == 26 and
            sum(r.expected_views) == |r.views|


    We also define a pseudo-set of bead pairs for each round
    - interactions \subsetof (Bead x Bead) : the set of interaction pairs processed so far in this round

    The interactions sub-set must only consist of interactions applied between elements in the round, and
    eventually it should sonsist of every possible pair of resident beads with resident and view beads.
    We will define the full interaction set as:
    - full_interactions(r) = { (b1,b2) : b1 \in r.resident, b2 \in ( r.resident + r.views ) }

    A round is interaction_complete if the following is true:
        interaction_complete(r) = r.interactions == full_interactions(r)
    This is currently defined in a way that allows interaction_complete to change over time,
    which is intentional. Essentially, if interaction_complete is false then there is compute
    work that could be done, while if it is true then either:
    - !recv_complete(r) -> We have no more compute work to do until another message arrives for this round; or
    - recv_complete(r) -> All compute work in the round is done (though there may still be sends to do if !send_complete(r))

    We are round_stable if there are no more messages to send or receive:
        round_stable(r) = send_complete(r) and recv_complete(r) and interaction_complete(r)

    Once a round is stable then:
    - That round will never change
    - If this is not an output round:
        - curr.outputs = { advance(b) for b in resident }
    - else:
        - curr.outputs = 0

    A round is complete if:
        round_complete(r) = round_stable(r) and r.outputs==0
    At this point we can move all the beads in 


    Now we break this into concrete messages and events (there are multiple ways this could be done).

    For sending we have an explicit ordering:
    - Migrate first : we want to get beads to their destinations, as the later they arrive,
      the more interactions have to be processed.
    - Share next : this is the bulk of the work
    - Expected : this is a message that only carries expected information, and is only
      needed if there were no migrations and no shares. So only needed if the cell is
      completely empty - this does happen, but is very rare.
    (Note that messages may arrive in a different order in the abstract hardware model)
    There is no particular ordering w.r.t. output messages.

    We send expected_migrate counts in the first message sent (regardless of type). We
    always know exactly how many beads are migrating at that point. If the first sent message
    is not of type migrate, then the migrate count must have been zero.

    We send expected_views in the first message sent after the cell becomes resident_complete.
    This may be any message type due to the ordering of all the messages coming in from
    neighbours, though this chain only extends one hop. At the start of each round every
    cell knows how many migration messages it produces, and that only needs to reach
    all neighbours (so no transitive dependency chains within a round).

    For interactions we enforce the invariant:
        r.interactions = { (a,b) : a in r.resident, b in (r.resident+r.viewed) }
    So any time we add a bead to the resident or the viewed set, we perform the
    interactions there and then. Recall that the resident and viewed sets are
    strictly disjoint, and that we expect all pair-wise interactions (including
    self-interactions) to be present. 

    attach_expected(r, msg):
        - if not r.sent_expected_migrations:
            - msg.expected_migrations = 0  // We could only have got here without sending expected migrations if there were zero migrations
        - if not r.sent_expected_views and resident_complete(r):
            - msg.expected_views = |r.resident|
            - r.sent_expected_views = 1

    process_expected(r, msg):
        - if msg.expected_migrations != None:
            - r.expected_migrations.add( msg.expected_migrations )
        - if msg.expected_views != None:
            - r.expected_views.add( msg.expected_views )

    Send/migrate: curr.migrating!=0
        - msg.type = migrate
        - msg.bead = pick(curr.migrating)
        - curr.migrating -= msg.bead
        - if not curr.sent_expected_migrations:
            - msg.expected_migrations = 1 + |curr.migrating|    // This must be the very first send
            - curr.sent_expected_migrations = true
        - attach_expected(r, msg)

    Send/share: curr.migrating==0 and curr.shared!=curr.resident  // We always send migrate before share
        - msg.type = share
        - msg.bead = pick( curr.resident - curr.shared )
        - curr.shared += msg.bead
        - attach_expected(r, msg)
    
    Send/expected: curr.migrating==0 and curr.shared==curr.resident
                    and ( !curr.sent_expected_migrations or (resident_complete(curr) && !curr.sent_expected_views)
        - msg.type = expected
        - msg.bead = none
        - attach_expected(r, msg)

    add_resident(r,a):
        - Pre: r.interactions = { (a,b) : a in r.resident, b in (r.resident+r.viewed) }
        - r.resident += a  // We include self-interaction
        - if is_output_round(r):
            - r.outputs += a
        - for b in r.viewed:
            - r.interactions += (a, b)
        - for b in r.resident:    
            - r.interactions += (a, b)
        -  Post: r.interactions = { (a,b) : a in r.resident, b in (r.resident+r.viewed) }

    Recv/migrate:
        - assert( t <= msg.t <= t+1)
        - dst = msg.t==t ? curr : next
        - if msg.bead is in cell:
            - add_resident( dst, msg.bead )
        - dst.migrations_received += 1
        - process_expected(dst, r)
    
    Recv/share:
        - assert( t <= msg.t <= t+1)
        - dst = msg.t==t ? curr : next
        - Pre: dst.interactions = { (a,b) : a in dst.resident, b in (dst.resident+dst.viewed) }
        - for a in dst.resident:
            - dst.interactions += (a, msg.bead)
        - dst.viewed += msg.bead
        - Post: dst.interactions = { (a,b) : a in dst.resident, b in (dst.resident+dst.viewed) }
        - process_expected(dst, r)

    Recv/expected:
        - assert( t <= msg.t <= t+1)
        - dst = msg.t==t ? curr : next
        - process_expected(dst, r)

    Step:  round_complete(r.curr) && t < T
        - Pre: next.interactions = { (a,b) : a in next.resident, b in (next.resident+next.viewed) }
        - for a in curr.resident:
            - a' = step(a)
            - if a' is in cell:
                add_resident( next, a' )
            - else:
                - next.migrating += a'
        - t = t+1
        
*/
struct GALSDPDEngineV1Handlers
{
    static constexpr size_t LOG2_MAX_BEADS_PER_CELL = 4;
    static constexpr size_t MAX_BEADS_PER_CELL = 1u<<LOG2_MAX_BEADS_PER_CELL;
    static constexpr size_t MAX_BEADS_PER_CELL_MASK = MAX_BEADS_PER_CELL-1;

    static constexpr size_t MAX_BEAD_TYPES=8;

    enum OutputFlags
    {
        RTS_INDEX_nhood=0,
        RTS_INDEX_output=1,

        RTS_FLAG_nhood=1<<RTS_INDEX_nhood,
        RTS_FLAG_output=1<<RTS_INDEX_output
    };

    struct bead_view_t
    {
        uint32_t id;
        float x[3];
        float v[3];
        int32_t src[3]; // Debug only
        uint32_t t; // Debug only
    };
    //static_assert(sizeof(bead_view_t)==28);
    
    struct bead_resident_t
        : bead_view_t
    {
        // Extras
        float f[3];
    };
    //static_assert(sizeof(bead_resident_t)==40);

    enum MessageType : uint8_t
    {
        None     = 0, // No message. Purely to do receive counts
        Migrate = 1,
        Share    = 2,  // Sharing view
        Output   = 3
    };

    enum EdgeMask : uint32_t
    {
        XLeft=1,
        XRight=2,
        YLeft=4,
        YRight=8,
        ZLeft=16,
        ZRight=32
    };

    struct message_t
    {
        uint32_t t;
        
        union{
            struct{
                MessageType type; // also uint8_t
                int8_t expected_migrations_count;        // If expected_count>=0, this is the number of messages sent in round t
                int8_t expected_shares_count;
                uint8_t _pad_;
            };
            uint32_t type_expected_raw;
        };
        
        bead_resident_t bead;
    };
    //static_assert(sizeof(message_t)<=56);

    struct round_state_t
    {
        uint32_t t;                    // time of this round

        uint8_t have_sent_expected_views;
        uint8_t have_sent_expected_migrations;

        uint32_t expected_views_received;
        uint32_t expected_views_delta;  // sum(expected_received) - num_views_received

        uint32_t expected_migrations_received;
        int32_t expected_migrations_delta; // sum(expected_migrations) - num_migrations_received

        uint32_t shared_count;   //  Beads [0,shared_count) are the shared set, beads [shared_count,num_resident) still need to be shared

        uint32_t num_resident;   //  [0,num_resident) are the resident set
        uint32_t num_migrating;   //  [MAX_BEADS_PER_CELL-num_migrating,MAX_BEADS_PER_CELL) are the migrating set
        
        // This is the viewed set
        uint32_t num_views;

        // These do not need to be cleared at start of round
        bead_resident_t beads[MAX_BEADS_PER_CELL];
        bead_view_t views[MAX_BEADS_PER_CELL*27];
    };

    static void clear_round_header(round_state_t &r)
    {
        uint32_t *p=(uint32_t*)&r;
        const size_t len=offsetof(round_state_t, beads);
        static_assert((len%4)==0);
        memzero32(p, len/4);
    }

    struct device_state_t
    {
        int32_t box[3];
        float dt;
        float scaled_inv_root_dt;
        struct {
            float conservative;
            float sqrt_dissipative;
        }interactions[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
        uint32_t max_t;

        uint32_t t;
        uint64_t t_hash;
        uint64_t t_seed;
        uint32_t interval_size;

        int32_t location[3];
        EdgeMask edge_mask;

        uint32_t rts;
 

        uint32_t intervals_todo;
        uint32_t interval_offset; // When offset reaches zero we output
        uint16_t outputs_waiting; // Number of outputs currently in the output queue
        uint16_t outputs_read_offset;   // Index of the first output in the queue
        
        round_state_t rounds[2];

        bead_resident_t outputs[MAX_BEADS_PER_CELL];
    };

    static bool is_resident_complete(const device_state_t &c, const round_state_t &r)
    {
        return c.t==r.t && r.expected_migrations_delta==0 && (r.expected_migrations_received==26);
    }

    static bool is_send_complete(const device_state_t &c, const round_state_t &r)
    {
        assert( r.have_sent_expected_views ? is_resident_complete(c,r) : true);
        assert( r.have_sent_expected_views ? r.have_sent_expected_migrations : true);
        return r.have_sent_expected_views && r.shared_count==r.num_resident;
    }

    static bool is_recv_complete(const device_state_t &c, const round_state_t &r)
    {
        return is_resident_complete(c,r) && r.expected_views_delta==0 && (r.expected_views_received==26);
    }

    static bool is_round_complete(const device_state_t &c, const round_state_t &r)
    {
        return is_send_complete(c,r) && is_recv_complete(c,r);
    }

    static bool is_output_queue_ready_for_push_burst(const device_state_t &c, unsigned num_beads)
    {
        return c.outputs_waiting + num_beads <= MAX_BEADS_PER_CELL;
    }

    static void output_queue_push(device_state_t &c, const bead_resident_t &bead)
    {
        assert(c.outputs_waiting < MAX_BEADS_PER_CELL);
        unsigned write_index=(c.outputs_read_offset+c.outputs_waiting) & MAX_BEADS_PER_CELL_MASK;
        memcpy32( c.outputs[write_index], bead );
        c.outputs_waiting += 1;
    }

    static void output_queue_pop(device_state_t &c, bead_resident_t &bead)
    {
        assert(c.outputs_waiting > 0);
        memcpy32(bead, c.outputs[c.outputs_read_offset]);
        c.outputs_read_offset = (c.outputs_read_offset+1) & MAX_BEADS_PER_CELL_MASK;
        c.outputs_waiting -=1;
    }

    
    static bool contains(const device_state_t &cell, const bead_resident_t &bead)
    {
        for(int i=0; i<3; i++){
            if(floor_nn(bead.x[i])!=cell.location[i]){
                return false;
            }
        }
        return true;
    }

    static bool advance_bead(const device_state_t &cell, bead_resident_t &b)
    {
        dpd_maths_core_half_step_raw::update_mom(cell.dt, b);
        dpd_maths_core_half_step_raw::update_pos_no_wrap(cell.dt, b);
        uint32_t em=cell.edge_mask;
        if(em){
            for(unsigned i=0; i<3; i++){
                if(em&0x3){
                    if( (em&0x1)){
                        if(b.x[i]<0){
                            b.x[i] += cell.box[i];
                        }
                    }else if(b.x[i]>=cell.box[i]){
                        b.x[i] -= cell.box[i];
                    }
                }
                em=em>>2;
            }
        }
        return contains(cell, b);
    }

    static void interact(const device_state_t &cell, uint32_t t, bead_resident_t &a, const bead_view_t &b, float *f_acc)
    {
        float dx[3];
        float dr2=0;
        for(int i=0; i<3; i++){
            dx[i] = a.x[i] - b.x[i];
            dr2 += dx[i]*dx[i];
            if(dr2 >= 1){
                return;
            }
        }

        /*fprintf(stderr, " %u : [%u,%u,%u]  %u <- %u,  [%f,%f,%f] <- [%f,%f,%f], p=%p\n",
            t, cell.location[0], cell.location[1], cell.location[2], a.id, b.id, a.x[0], a.x[1], a.x[2], b.x[0], b.x[1], b.x[2], f_acc
        );*/

        auto inter=cell.interactions[ BeadHash{a.id}.get_bead_type()*MAX_BEAD_TYPES + BeadHash{b.id}.get_bead_type() ];

        float f[3];
        dpd_maths_core_half_step_raw::calc_force<false,float,float[3],float[3]>(
            cell.scaled_inv_root_dt,
            cell.t_hash,
            dx,
            pow_half(dr2),
            0.0f,
            0.0,
            inter.conservative,
            inter.sqrt_dissipative,
            BeadHash{a.id},
            BeadHash{b.id},
            a.v,
            b.v,
            f
        );

        vec3_add(a.f, f);
        if(f_acc){
            vec3_sub(f_acc, f);
        }
    }

    static void add_resident(const device_state_t &cell, round_state_t &r, const bead_resident_t &b)
    {
        for(unsigned i=0; i<r.num_views; i++){
            assert(b.id != r.views[i].id);
        }
        for(unsigned i=0; i<r.num_resident; i++){
            assert(b.id != r.beads[i].id);
        }
        for(unsigned i=0; i<r.num_migrating; i++){
            assert(b.id != r.beads[MAX_BEADS_PER_CELL-1-i].id);
        }

        assert(b.t==r.t);

        auto &rb=r.beads[r.num_resident];
        memcpy32(rb, b);
        float *f=rb.f;
        assert(vec3_is_zero(f));
        for(unsigned i=0; i<r.num_resident; i++){
            interact( cell, r.t, r.beads[i], rb, f);
        }
        for(unsigned i=0; i<r.num_views; i++){
            interact( cell, r.t, rb, r.views[i], nullptr);
        }
        r.num_resident += 1;
    }

    static void add_migrating(const device_state_t &c, round_state_t &r, const bead_resident_t &b)
    {
        assert(b.t==r.t);

        auto &rb=r.beads[MAX_BEADS_PER_CELL-1-r.num_migrating];
        memcpy32(rb, b);
        assert(vec3_is_zero(rb.f));
        r.num_migrating += 1;
    }

    static void add_view(const device_state_t &cell, round_state_t &r, const bead_resident_t &b)
    {


        /*fprintf(stderr, " %u : [%u,%u,%u]  add_view(%u) from [%u,%u,%u], num_views=%u, dp=%p, rp=%p\n",
            r.t, cell.location[0], cell.location[1], cell.location[2], b.id,
            b.src[0], b.src[1], b.src[2],
            r.num_views,
            &cell, &r
        );*/

        assert(b.t==r.t);

        auto &rv=r.views[r.num_views];
        memcpy32(rv, (bead_view_t&)b);

        for(unsigned i=0; i<r.num_views; i++){
            assert(rv.id != r.views[i].id);
        }
        for(unsigned i=0; i<r.num_resident; i++){
            assert(rv.id != r.beads[i].id);
        }
        for(unsigned i=0; i<r.num_migrating; i++){
            assert(rv.id != r.beads[MAX_BEADS_PER_CELL-1-i].id);
        }

        // If needed, apply wrap-around.
        // optimised for case where no wrapping needed
        uint32_t em=cell.edge_mask;
        if(em){
            // Vast majority of cells never take this path
            for(unsigned i=0; i<3; i++){
                if(em &0x3){
                    if( em&0x1 ){
                        if(rv.x[i]>2){
                            rv.x[i] -= cell.box[i];
                        }
                    }else if(rv.x[i]<=2){
                        assert(em&0x2);
                        rv.x[i] += cell.box[i];
                    }
                }
                em>>=2;
            }
        }

        for(unsigned i=0; i<r.num_resident; i++){
            interact(cell, r.t, r.beads[i], rv, nullptr);
        }
    }

    static bool is_output_round(const device_state_t &c, uint32_t t)
    {
        assert(c.t <= t && t<=c.t+1 );
        // This function is not on a very hot path (only called once
        // per round), but it could be more efficient
        if(c.t==t){
            return c.interval_offset==0;
        }
        assert(c.t+1==t);
        return (c.interval_offset==1 && c.interval_size==1);
    }

    static void round_invariants(const device_state_t &cell, const round_state_t &r)
    {
        if(is_round_complete(cell,r)){
            assert( is_round_complete(cell,r) );
            assert( is_send_complete(cell,r) );
            assert( is_recv_complete(cell,r) );
        }

        for(unsigned i=0; i<r.num_resident; i++){
            assert(r.beads[i].t==r.t);
            assert(contains(cell, r.beads[i]));

            for(unsigned j=i+1; j<r.num_resident; j++){
                assert(r.beads[i].id != r.beads[j].id);
            }
            for(unsigned j=0; j<r.num_migrating; j++){
                assert(r.beads[i].id != r.beads[MAX_BEADS_PER_CELL-1-j].id);
            }
            for(unsigned j=0; j<r.num_views; j++){
                assert(r.beads[i].id != r.views[j].id);
            }
        }
    }

    static void invariants(const device_state_t &cell)
    {
        round_invariants(cell, cell.rounds[0]);
        round_invariants(cell, cell.rounds[1]);
        assert(cell.t <= cell.max_t);

        const auto &curr=cell.rounds[cell.t&1];
        assert(curr.t==cell.t);

        if(cell.t <= cell.max_t){
            assert_implies(curr.have_sent_expected_views, curr.have_sent_expected_migrations);
            assert_implies(curr.num_migrating, cell.rts == RTS_FLAG_nhood);
            assert_implies(curr.shared_count < curr.num_resident, cell.rts == RTS_FLAG_nhood );
            assert_implies(!curr.have_sent_expected_migrations, cell.rts == RTS_FLAG_nhood );
            assert_implies(!curr.have_sent_expected_views && is_resident_complete(cell,curr), cell.rts == RTS_FLAG_nhood);
        }else{
            assert(cell.t+1==cell.max_t);
            assert(curr.have_sent_expected_migrations==0);
            assert(curr.have_sent_expected_views==0);
            assert(curr.shared_count==0);
        }
        
        if(cell.rts!=RTS_FLAG_nhood){
            assert_implies(cell.outputs_waiting>0, cell.rts == RTS_FLAG_output);
        }
    }

    static void advance(device_state_t &c)
    {
        auto &curr=c.rounds[c.t&1];
        auto &next=c.rounds[1-(c.t&1)];

        assert(is_round_complete(c, curr));
        assert(curr.t==c.t);
        assert(next.t==c.t+1);

        c.t += 1;

        bool add_to_output_queue=is_output_round(c, c.t);

        //fprintf(stderr, "  advance : %u -> %u, is_output=%d, resident=%u\n", c.t-1, c.t, add_to_output_queue, curr.num_resident);

        for(unsigned i=0; i<curr.num_resident; i++){
            auto &b=curr.beads[i];
            assert(b.t==curr.t);
            
            b.t += 1;
            if(add_to_output_queue){
                output_queue_push(c, b);
            }

            /*fprintf(stderr, " %u -> %u : x=[%f,%f,%f], v=[%f,%f,%f], f=[%f,%f,%f]\n",
                c.t, c.t+1, b.x[0], b.x[1], b.x[2], b.v[0], b.v[1], b.v[2],
                b.f[0], b.f[1], b.f[2]
            );*/

            bool still_resident=advance_bead(c, b);

            if(still_resident){
                add_resident(c, next, b);
            }else{
                add_migrating(c, next, b);
            }
        }

        clear_round_header(curr);
        curr.t=c.t+1; // curr is now "next"

        if(c.interval_offset==0){
            c.interval_offset=c.interval_size;
        }
        c.interval_offset -= 1;

        if(c.t > c.max_t){
            // We have halted
            assert(add_to_output_queue); // Why no output in final time step?
            if(c.outputs_waiting>0){
                c.rts=RTS_FLAG_output;
            }else{
                c.rts=0;
            }
        }else{
            // It is guaranteed we are ready to send, regardless of
            // the number of resident
            c.rts = RTS_FLAG_nhood;
        }

        /*fprintf(stderr, " %u : [%u,%u,%u], advanced, new_views=%u, new_resident=%u\n",
            c.t, c.location[0], c.location[1], c.location[2], next.num_views, next.num_resident
        );*/ 
    }


    static void update_rts_and_maybe_step(device_state_t &cell, const round_state_t &curr)
    {
        assert(&curr==&cell.rounds[cell.t&1]);
        if(curr.num_migrating
            || curr.shared_count!=curr.num_resident  // fast paths
            || !curr.have_sent_expected_migrations
            || ( is_resident_complete(cell,curr) && !curr.have_sent_expected_views) )
        {
            cell.rts = RTS_FLAG_nhood;

            invariants(cell);
            return;
        }
        
        cell.rts = 0;

        if(curr.have_sent_expected_views){
            // This could mean we are ready to step
            assert(curr.num_migrating==0);
            assert(curr.shared_count==curr.num_resident);
            assert(curr.have_sent_expected_migrations);
            assert(curr.have_sent_expected_views);

            if(curr.expected_views_received==26 && curr.expected_views_delta==0
                && curr.expected_migrations_received==26 && curr.expected_migrations_delta==0){
                // We are round complete
                assert(is_round_complete(cell,curr));
                if(cell.t<cell.max_t && is_output_queue_ready_for_push_burst(cell, curr.num_resident)){
                    advance(cell);
                    invariants(cell);
                    return;
                }
            }
        }

        if(cell.outputs_waiting){
            cell.rts=RTS_FLAG_output;
        }

        invariants(cell);
    }

    static uint32_t calc_rts(const device_state_t &cell)
    {
        return cell.rts;
    }

    static void on_init(device_state_t &cell)
    {
        invariants(cell);
    }

    static bool on_barrier(device_state_t &cell)
    {
        return true;
    }

    static void on_send(device_state_t &cell, message_t &msg)
    {
        invariants(cell);

        assert(cell.rts);
        round_state_t &curr=cell.rounds[cell.t&1];

        // Debug only
        for(int i=0; i<3; i++){
            msg.bead.src[i]=cell.location[i];
        }

        if(cell.rts==RTS_FLAG_nhood){
            msg.t=curr.t;
            msg.expected_shares_count=-1;
            msg.expected_migrations_count=-1;

            if(curr.num_migrating){
                assert(cell.rts==RTS_FLAG_nhood);
                msg.type=Migrate;
                memcpy32(msg.bead, curr.beads[MAX_BEADS_PER_CELL-curr.num_migrating]);
                if(!curr.have_sent_expected_migrations){
                    msg.expected_migrations_count = curr.num_migrating;
                    curr.have_sent_expected_migrations = 1;
                }
                curr.num_migrating -= 1;
            }else if(curr.shared_count != curr.num_resident){
                assert(cell.rts==RTS_FLAG_nhood);
                msg.type=Share;
                memcpy32(msg.bead, curr.beads[curr.shared_count]);
                curr.shared_count += 1;
            }else if(!curr.have_sent_expected_migrations){
                msg.type=None;
            }else if(!curr.have_sent_expected_views && is_resident_complete(cell,curr)){
                msg.type=None;
            }else{
                assert(0);
            }

            if(!curr.have_sent_expected_migrations){
                msg.expected_migrations_count=0;
                curr.have_sent_expected_migrations=1;
            }
            if(!curr.have_sent_expected_views && is_resident_complete(cell,curr)){
                msg.expected_shares_count=curr.num_resident;
                curr.have_sent_expected_views=1;
            }

            update_rts_and_maybe_step(cell, curr);

            invariants(cell);
        }else{
            assert(cell.rts==RTS_FLAG_output);
            assert(cell.outputs_waiting);
            msg.type=Output;
            output_queue_pop(cell, msg.bead);

            msg.t=msg.bead.t;

            if(cell.t < cell.max_t){
                // Still running normallu
                update_rts_and_maybe_step(cell, curr);
                invariants(cell);
            }else if(cell.outputs_waiting==0){
                // Finished and just drained final output. Go completely quiet
                cell.rts=0;
                invariants(cell);
            }
        }

        invariants(cell);
    }

    static void on_recv(device_state_t &cell, const message_t &msg)
    {
        invariants(cell);

        assert(cell.t <= msg.t && msg.t <= cell.t+1);
        auto &dst = cell.rounds[msg.t&1];
        assert(dst.t==msg.t);

        if(msg.type==Share){
            add_view(cell, dst, msg.bead);
            dst.expected_views_delta -= 1;
        }else if(msg.type==Migrate){
            assert(!is_resident_complete(cell,dst));
            if(contains(cell, msg.bead)){
                add_resident(cell, dst, msg.bead);
            }
            dst.expected_migrations_delta -= 1;
        }else{
            assert(msg.type==None);
        }

        if(msg.expected_migrations_count>=0){
            dst.expected_migrations_received += 1;
            dst.expected_migrations_delta += msg.expected_migrations_count;
        }
        if(msg.expected_shares_count>=0){
            dst.expected_views_received +=1;
            dst.expected_views_delta += msg.expected_shares_count;
        }

        if(msg.t==cell.t){
            update_rts_and_maybe_step(cell, dst);
        }

        invariants(cell);
    }


};

#endif
