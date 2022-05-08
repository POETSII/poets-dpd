#ifndef basic_dpd_engine_v8_raw_handlers_hpp
#define basic_dpd_engine_v8_raw_handlers_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw_handlers.hpp"

struct BasicDPDEngineV8RawHandlers
    : public BasicDPDEngineV5RawHandlers
{
    struct raw_bead_share_t
    {
        raw_bead_view_t beads[2];
    };
    static_assert(sizeof(raw_bead_share_t) <= 60);

    struct device_state_t
        : BasicDPDEngineV5RawHandlers::device_state_t
    {
        uint32_t loc;
    };


    template<bool EnableLogging>
    static __attribute__((noinline)) bool on_barrier(device_state_t &cell)
    {
        switch(cell.phase){
            default: assert(false); // fallthrough
            case PreMigrate:
            case Outputting: 
            case SharingAndForcing: return on_barrier_pre_migrate(cell); break;
            case Migrating: return on_barrier_pre_share<EnableLogging>(cell); break;            
        }
    }

    static void /*__attribute__ ((noinline)) __attribute__((optimize("O0")))*/ on_init(device_state_t &cell)
    {
        BasicDPDEngineV5RawHandlers::on_init(cell);
    }

    template<bool EnableLogging>
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

    static void on_send_share(device_state_t &cell, raw_bead_share_t &outgoing)
    {
        auto resident=make_bag_wrapper(cell.resident);

        assert(cell.share_todo>0);

        raw_bead_view_t *begin=outgoing.beads;
        raw_bead_view_t *end=outgoing.beads+std::size(outgoing.beads);

        while(cell.share_todo && begin!=end){
            --cell.share_todo;
            const auto &b = resident[cell.share_todo];
            copy_bead_view::copy( begin, &b );
            begin++;
            break;
        }
        while(begin!=end){
            begin->id=0xFFFFFFFFul;
            begin++;
        }

        if(cell.share_todo==0){
            calc_intra_forces<false>(cell);

            cell.rts &= ~RTS_FLAG_share;
        }
    }

    enum BondLevel
    {
        NotBonded=0,
        SpatialBonding,
        HookeanBonding // also implies spatial bonding
    };

    template<bool EnableLogging>
    static int  interact(
        device_state_t &cell,
        raw_bead_resident_t &bead, const float *bead_x, uint32_t incoming_id, const float *incoming_x, const float *incoming_v,
        float *f3
        )
    {
        float dx[3];
        vec3_sub(dx, bead.x, incoming_x);
        float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        if(dr_sqr >=1 || dr_sqr < MIN_DISTANCE_CUTOFF_SQR){ // The min threshold avoid large forces, and also skips self-interaction
            return NotBonded;
        }

        assert(bead.id != incoming_id);

        float dr=pow_half(dr_sqr);

        BondLevel bondLevel=SpatialBonding;
        float kappa=0.0f;
        float r0=cell.bond_r0;
        if(!(is_monomer(incoming_id) || is_monomer(bead.id))){
            if(get_polymer_id(bead.id) == get_polymer_id(incoming_id) ){
                auto other_polymer_offset=get_polymer_offset(incoming_id);
                for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
                    if(other_polymer_offset==bead.bond_partners[i]){
                        bondLevel=HookeanBonding;
                        kappa=cell.bond_kappa;
                        break;
                    }
                }
            }
        }

        auto bead_type1=get_bead_type(bead.id);
        auto bead_type2=get_bead_type(incoming_id);
        auto interactions=cell.interactions[MAX_BEAD_TYPES*bead_type1+bead_type2];

        float incoming_vv[3]={incoming_v[0],incoming_v[1],incoming_v[2]};
        dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float*>(
            cell.scaled_inv_root_dt,
            cell.t_hash,
            dx, dr,
            kappa, r0, 
            interactions.conservative,
            interactions.sqrt_dissipative,
            get_hash_code(bead.id), get_hash_code(incoming_id),
            bead.v, incoming_vv,
            f3
        );

        vec3_add(bead.f, f3);

        return bondLevel;
    }

    static void __attribute__((noinline)) cache_bond(device_state_t &cell, unsigned bead_i, raw_bead_resident_t &bead, uint32_t incoming_id, const float *neighbour_x)
    {
        static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
        bool hit = get_polymer_offset(incoming_id) == bead.angle_bonds[0].partner_head || get_polymer_offset(incoming_id) == bead.angle_bonds[0].partner_tail;
        if(!hit){
            return;
        }

        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);
        
        //std::cerr<<"K: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
        unsigned cached_bond_index=cached_bonds.size();
        raw_cached_bond_t tmp;
        tmp.bead_hash=BeadHash{incoming_id}.hash;
        vec3_copy(tmp.x, neighbour_x);
        cached_bonds.push_back(tmp);
        //std::cerr<<" caching, bead_i="<<bead_i<<", bead_id="<<get_hash_code(bead.id)<<" target="<<get_hash_code(incoming.id)<<"\n";

        //std::cerr<<"  CB: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
        if(cell.cached_bond_indices[bead_i]==0xFF){
            //std::cerr<<"  S\n";
            cell.cached_bond_indices[bead_i] = cached_bond_index;
        }else{
            // Once both partners have arrived, we calculate force and push onto outgoing_forces
            //std::cerr<<"Force!`\n";
            assert(bead.id == make_bag_wrapper(cell.resident)[bead_i].id);
            calc_angle_force_for_middle_bead(cell, bead, cell.cached_bond_indices[bead_i], cached_bond_index);
            cell.rts |= RTS_FLAG_force;
            //std::cerr<<"DoenFr\n";
        }
    }

    template<bool EnableLogging>
    static void on_recv_share(device_state_t &cell, const raw_bead_share_t &incoming_share)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<")\n";

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);

        auto begin=incoming_share.beads;
        auto end=incoming_share.beads+std::size(incoming_share.beads);

        while(begin!=end){
            const auto &incoming=*begin;
            if(incoming.id==0xFFFFFFFFul){
                break;
            }
            ++begin;
            
            float neighbour_x[3];
            vec3_copy(neighbour_x, incoming.x);
            int32_t neighbour_cell_pos[3];
            vec3_floor_nn(neighbour_cell_pos, neighbour_x);
            
            // TODO : don't send these messages to self
            if(neighbour_cell_pos[0]==cell.location[0] && neighbour_cell_pos[1]==cell.location[1] && neighbour_cell_pos[2]==cell.location[2]){
                //std::cerr<<" Skipping self.\n";
                return;
            }

            for(unsigned i=0; i<resident.size(); i++){
                assert(incoming.id != resident[i].id);
            }
            

            for(int d=0; d<3; d++){
                if(cell.location[d]==0 && neighbour_cell_pos[d]==cell.box[d]-1){
                    neighbour_x[d] -= cell.box[d];
                }else if(cell.location[d]==cell.box[d]-1 && neighbour_cell_pos[d]==0){
                    neighbour_x[d] += cell.box[d];
                }
            }

            unsigned cached_bond_index=cached_bonds.size(); // This is the index it will have, _if_ it is cached
            for(unsigned bead_i=0; bead_i < resident.size(); bead_i++){
                auto &bead=resident[bead_i];

                float f[3];
                auto bondLevel=interact<EnableLogging>(cell, bead, bead.x, incoming.id, neighbour_x, incoming.v, f );
                if(bondLevel==HookeanBonding){
                    cache_bond(cell, bead_i, bead, incoming.id, neighbour_x);
                }
            }
        }
    }

    template<bool EnableLogging>
    static void calc_intra_forces(device_state_t &cell)
    {
        auto resident=make_bag_wrapper(cell.resident);

        for(unsigned i=0; i+1<resident.size(); i++){
            for(unsigned j=1; j<resident.size(); j++){
                auto &bead=resident[i];
                auto &neighbour=resident[j];

                float f[3];
                auto bondLevel=interact<EnableLogging>(cell, bead, bead.x, neighbour.id, neighbour.x, neighbour.v, f);
                //std::cerr<<"bondLevel="<<bondLevel<<", f="<<f[0]<<","<<f[1]<<","<<f[2]<<"\n";
                if(bondLevel > NotBonded){
                    vec3_sub(neighbour.f, f);
                }

                if(bondLevel == HookeanBonding){
                    cache_bond(cell, i, bead, neighbour.id, neighbour.x);
                    cache_bond(cell, j, neighbour, bead.id, bead.x);
                }
            }
        }
    }

};

#endif
