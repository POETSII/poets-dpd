

/*
While writing this it became unclear why it is needed for the
set of polymers describeable by the PolymerType as it stands.
BondPairs are explicitly represented in terms of Bonds, which
means recalculation can be avoided quite easily.

There also seem to be certain polymers that can't be described
by the current PolymerType, for example:

    D 
   /|
A-B-C

If we represent it as:
 A-B-C as (A->B)->(B->C)
 A-B-D as (A->B)->(B->D)
 B-C-D as (B->C)->(C->D)
then we have a problem for B-D-C
 (B->D)->(D->C)
 or
 (C->D)->(D->B)
either way we end up with reversed bonds on one leg.
This could be handled by changing theta0 I guess, but
not sure whether it really works.
*/
#error "Not tested, also not clear if needed."


/*!
The purpose of this class is to take a PolymerType and turn
it into something that can be efficiently evaluated at run-time.
Exactly which optimisations it performs will probably change
over time, but it should at least reduce repeated evaluations
of bond distances and lengths.

*/
struct bond_jit
{
    struct BeadValue
    {
        intptr_t handle;
        BeadHash hash;
        vec3f_t x;
        vec3f_t f;
    };


    struct jit_plan
    {
        struct PairInfo
        {
            uint16_t offset0;
            uint16_t offset1;
            float r0;
            float kappa;
        };
        struct TripleInfo
        {
            uint16_t offset0;
            int16_t offset1;   // If negative then we need to do -dx
            // Do not support non-straight bonds yet
            float kappa;
            uint16_t dst_offset0; // Destinations for the three forces
            uint16_t dst_offset1;
            uint16_t dst_offset2;
        };

        unsigned bead_count;
        std::vector<PairInfo> pair_info;
        std::vector<TripleInfo> triple_info;
    };

    jit_plan build_jit_plan(const PolymerType &pt)
    {
        jit_plan res;
        res.bead_count=pt.bead_types.size();

        // Map to positive number if it is directly there, negative if it is flipped
        std::map<std::pair<unsigned,unsigned>,int> bead_pair_locations;
        for(const auto &bp : pt.bonds){
            auto it=bead_pair_locations.find({bp.bead_offset_head,bp.bead_offset_tail});
            if(it==bead_pair_locations.end()){
                bead_pair_locations.insert( {{bp.bead_offset_head,bp.bead_offset_tail},res.pair_info.size() } );
                bead_pair_locations.insert( {{bp.bead_offset_tail,bp.bead_offset_head},-(int)res.pair_info.size() } );
                res.bead_pairs.push_back({bp.bead_offset_head, bp.bead_offset_tail, bp.r0, bp.kappa });
            }
        }

        for(const auto &bp : pt.bond_pairs){
            const auto &bh=pt.bonds[bp.bond_pair_head];
            const auto &bt=pt.bonds[bp.bond_pair_tail];
            unsigned o0=bh.bead_offset_head, o1=bh.bead_offset_tail, o2=bt.bead_offset_tail;
            int l0=bead_pair_locations.find({o0,o1});
            int l1=bead_pair_locations.find({o1,o2});
            if(l0<0){
                l0=-l0;
                l1=-l1;
            }  
        }
    }

    struct jit_context
    {
        struct PairValue
        {
            vec3r_t dx; 
            float r;
        };

        std::vector<BeadPos> bead_positions;
        std::vector<PairValue> bead_pairs;

        template<TFetchBead fetch_bead>
        void execute(
            const jit_plan &plan,
            uint32_t polymer_id,
            TReader &&reader,  // reader(unsigned polymer_id, unsigned polymer_offset) -> BeadValue
            TWriter &&writer   // writer(BeadValue b) -> void
        ){
            assert(bead_positions.size() >= plan.bead_hashes.size());
            assert(position_deltas.size() >= plan.position_delta_pairs.size());

            for(unsigned i=0; i<plan.bead_hashes.size(); i++){
                bead_positions[i] = reader( polymer_id, i );
            }
            for(unsigned i=0; i<plan.position_delta_pairs.size(); i++){
                const PairInfo &pi=plan.pair_info[i];
                PairValue val;
                
                val.dx = bead_positions[pi.offset0] - bead_positions[pi.offset1];
                val.r = val.dx.l2_norm();

                float hookean=(pi.r0-dr)*pi.kappa;
                vec3f_t force=dx * (hookean/val.r);

                bead_positions[pi.offset0].f += force;
                bead_positions[pi.offset1].f -= force;

                bead_pairs[i] = val; // TODO: What if there are no triples
            }
            for(unsigned i=0; i<plan.triple_info.size(); i++){
                const TripleInfo &ti=plan.triple_info[i];
                const auto &pair0 = bead_pairs[ti.offset0];
                const auto &pair1 = bead_pairs[ti.offset1];

                vec3r_t pair1_dx = ti.offset1 < 0 ? -pair1.dx : pair1.dx;

                float magProduct = pair0.r*pair1.r;
                if(magProduct > 0.0001f)
                {
                    float b1MagSq		= pair0.r*pair0.r;
                    float b2MagSq		= pair1.r*pair1.r;
                    float b1Dotb2		= pair0.dx[0]*pair1_dx[0] + pair0.dx[1]*pair1_dx[1] + pair0.dx[2]*pair1_dx[2];
                    float b1b2Overb1Sq	= b1Dotb2/b1MagSq;
                    float b1b2Overb2Sq	= b1Dotb2/b2MagSq;
                    float cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

                    float forceMag = ti.kappa/magProduct;

                    vec3f_t headForce =((pair0.dx*b1b2Overb1Sq)-pair1.dx) * forceMag;
                    vec3f_t tailForce = (pair0.dx- (pair1.dx*b1b2Overb2Sq)) * forceMag;
                    vec3f_t middleForce = -( headForce + tailForce );

                    bead_positions[ti.dst_offset0] += headForce;
                    bead_positions[ti.dst_offset1] += middleForce;
                    bead_positions[ti.dst_offset2] += tailForce;
                }
            }
            for(unsigned i=0; i<bead_positions.size(); i++){
                writer(bead_positions[i]);
            }
        }
    };


};