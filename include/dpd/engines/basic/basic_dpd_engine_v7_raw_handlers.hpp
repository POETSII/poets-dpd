#ifndef basic_dpd_engine_v7_raw_handlers_hpp
#define basic_dpd_engine_v7_raw_handlers_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw_handlers.hpp"

struct BasicDPDEngineV7RawHandlers
    : public BasicDPDEngineV5RawHandlers
{
    static const bool USE_X_CACHE=true;

    struct raw_bead_share_t
    {
        raw_bead_view_t beads[2];
    };
    static_assert(sizeof(raw_bead_share_t) <= 56);

    struct device_state_t
        : BasicDPDEngineV5RawHandlers::device_state_t
    {
        float x_cache[MAX_BEADS_PER_CELL*3];
    };


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

    static bool on_barrier_pre_share(device_state_t &cell)
    {
        assert(cell.phase==Migrating);

        auto resident=make_bag_wrapper(cell.resident);
        assert(make_bag_wrapper(cell.force_outgoing).empty());
        assert(make_bag_wrapper(cell.migrate_outgoing).empty());
        assert(cell.share_todo==0);

        cell.share_todo = resident.size();

        float *x_cache=cell.x_cache; 
        for(unsigned i=0; i<resident.size(); i++, x_cache+=3){ 
            cell.cached_bond_indices[i]=0xff;
            if(USE_X_CACHE){
                for(int j=0; j<3; j++){
                    x_cache[j] = resident[i].x[j];
                }
            }
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
            copy_bead_view( begin, &b );
            begin++;
        }
        while(begin!=end){
            begin->id=0xFFFFFFFFul;
            begin++;
        }

        if(cell.share_todo==0){
            cell.rts &= ~RTS_FLAG_share;
        }
    }

    static void on_recv_share(device_state_t &cell, const raw_bead_share_t &incoming_share)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<"), p="<<&cell<<", nres="<<cell.resident.n<<", other="<<get_hash_code(incoming.id)<<"\n";

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
            for(int d=0; d<3; d++){
                if(cell.location[d]==0 && neighbour_cell_pos[d]==cell.box[d]-1){
                    neighbour_x[d] -= cell.box[d];
                }else if(cell.location[d]==cell.box[d]-1 && neighbour_cell_pos[d]==0){
                    neighbour_x[d] += cell.box[d];
                }
            }

            bool cached=false;

            const float *bead_x=cell.x_cache;
            unsigned cached_bond_index=cached_bonds.size(); // This is the index it will have, _if_ it is cached
            for(unsigned bead_i=0; bead_i < resident.size(); bead_i++, bead_x+=3){
                auto &bead=resident[bead_i];

                //fprintf(stderr, "    Recv : bead_i=%u, home=%u, other=%u\n", bead_i, get_hash_code(bead.id), get_hash_code(incoming.id));

                // This implicitly interacts each bead with itself, which is handled with a
                // distance check in calc_force.
                float dx[3];
                if(USE_X_CACHE){
                    vec3_sub(dx, bead_x, neighbour_x);
                }else{
                    vec3_sub(dx, bead.x, neighbour_x);
                }
                float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr_sqr >=1 || dr_sqr < float(1e-5)){ // The min threshold avoid large forces, and also skips self-interaction
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
                float conStrength=cell.conservative[MAX_BEAD_TYPES*bead_type1+bead_type2];

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

                if(kappa!=0.0f){
                    //std::cerr<<"K: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
                    if(!cached){
                        raw_cached_bond_t tmp;
                        tmp.bead_hash=get_hash_code(incoming.id);
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
    }

};

#endif
