#ifndef basic_dpd_engine_v7_raw_handlers_hpp
#define basic_dpd_engine_v7_raw_handlers_hpp

#include "dpd/engines/basic/basic_dpd_engine_v5_raw_handlers.hpp"

template<bool EnableLogging=false, bool USE_X_CACHE=false>
struct BasicDPDEngineV7RawHandlers
    : public BasicDPDEngineV5RawHandlers
{
    using Base  = BasicDPDEngineV5RawHandlers;    

    struct raw_bead_share_t
    {
        raw_bead_view_t beads[2];
    };
    static_assert(sizeof(raw_bead_share_t) <= 56);

    struct device_state_t
        : Base::device_state_t
    {
        bool is_edge;
        float x_cache[MAX_BEADS_PER_CELL*3];
    };


    static bool on_barrier(device_state_t &cell)
    {
        switch(cell.phase){
            default: assert(false); // fallthrough
            case Phase::PreMigrate:
            case Phase::Outputting: 
            case Phase::SharingAndForcing: return Base::on_barrier_pre_migrate(cell); break;
            case Phase::Migrating: return on_barrier_pre_share(cell); break;            
        }
    }

    static bool on_barrier_pre_share(device_state_t &cell)
    {
        assert(cell.phase==Phase::Migrating);

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

        cell.phase=Phase::SharingAndForcing;
        cell.rts=cell.share_todo==0 ? 0 : OutputFlags::RTS_FLAG_share;

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
            Base::copy_bead_view::copy( begin, &b );
            begin++;
        }
        while(begin!=end){
            begin->id=0xFFFFFFFFul;
            begin++;
        }

        if(cell.share_todo==0){
            cell.rts &= ~OutputFlags::RTS_FLAG_share;
        }
    }

    static void on_recv_share(device_state_t &cell, const raw_bead_share_t &incoming_share)
    {
        //std::cerr<<"  Recv: ("<<cell.location[0]<<","<<cell.location[1]<<","<<cell.location[2]<<"), p="<<&cell<<", nres="<<cell.resident.n<<", other="<<get_hash_code(incoming.id)<<"\n";

        auto resident=make_bag_wrapper(cell.resident);
        auto cached_bonds=make_bag_wrapper(cell.cached_bonds);

        auto begin=incoming_share.beads;
        auto end=incoming_share.beads+std::size(incoming_share.beads);

        bool is_edge=true; //cell.is_edge;

        while(begin!=end){
            const auto &incoming=*begin;
            if(incoming.id==0xFFFFFFFFul){
                break;
            }
            ++begin;
            
            float neighbour_x[3];
            vec3_copy(neighbour_x, incoming.x);

            if(is_edge){
                int32_t neighbour_cell_pos[3];
                vec3_floor_nn(neighbour_cell_pos, neighbour_x);
                for(int d=0; d<3; d++){
                    if(cell.location[d]==0 && neighbour_cell_pos[d]==cell.box[d]-1){
                        neighbour_x[d] -= cell.box[d];
                    }else if(cell.location[d]==cell.box[d]-1 && neighbour_cell_pos[d]==0){
                        neighbour_x[d] += cell.box[d];
                    }
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
                    for(int d=0;d<3;d++){
                        assert(bead_x[d]==bead.x[d]);
                    }
                    vec3_sub(dx, bead_x, neighbour_x);
                }else{
                    vec3_sub(dx, bead.x, neighbour_x);
                }
                float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr_sqr >=1 || dr_sqr < MIN_DISTANCE_CUTOFF_SQR){ // The min threshold avoid large forces, and also skips self-interaction
                    continue;
                }
                if(bead.id==incoming.id){
                    continue;
                }

                float dr=pow_half(dr_sqr);

                float kappa=0.0f;
                float r0=cell.bond_r0;
                if(!(BeadHash{incoming.id}.is_monomer() || BeadHash{bead.id}.is_monomer())){
                    if(BeadHash{bead.id}.get_polymer_id() == BeadHash{incoming.id}.get_polymer_id() ){
                        auto other_polymer_offset=BeadHash{incoming.id}.get_polymer_offset();
                        for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
                            if(other_polymer_offset==bead.bond_partners[i]){
                                kappa=cell.bond_kappa;
                                break;
                            }
                        }
                    }
                }

                auto bead_type1=BeadHash{bead.id}.get_bead_type();
                auto bead_type2=BeadHash{incoming.id}.get_bead_type();
                auto strength=cell.interactions[MAX_BEAD_TYPES*bead_type1+bead_type2];

                float f[3];
                dpd_maths_core_half_step_raw::calc_force<EnableLogging,float,float[3],float[3]>(
                    cell.scaled_inv_root_dt,
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
                        tmp.bead_hash=incoming.id;
                        vec3_copy(tmp.x, neighbour_x);
                        cached_bonds.push_back(tmp);
                        //std::cerr<<" caching, bead_i="<<bead_i<<", bead_id="<<get_hash_code(bead.id)<<" target="<<get_hash_code(incoming.id)<<"\n";
                        cached=true;
                    }

                    static_assert(MAX_ANGLE_BONDS_PER_BEAD==1);
                    bool hit = BeadHash{incoming.id}.get_polymer_offset() == bead.angle_bonds[0].partner_head || BeadHash{incoming.id}.get_polymer_offset() == bead.angle_bonds[0].partner_tail;
                    if(hit){
                        //std::cerr<<"  CB: "<<get_hash_code(bead.id)<<" - "<<get_hash_code(incoming.id)<<"\n";
                        if(cell.cached_bond_indices[bead_i]==0xFF){
                            //std::cerr<<"  S\n";
                            cell.cached_bond_indices[bead_i] = cached_bond_index;
                        }else{
                            // Once both partners have arrived, we calculate force and push onto outgoing_forces
                            //std::cerr<<"Force!`\n";
                            assert(bead.id == resident[bead_i].id);
                            Base::calc_angle_force_for_middle_bead(cell, bead, cell.cached_bond_indices[bead_i], cached_bond_index);
                            cell.rts |= OutputFlags::RTS_FLAG_force;
                            //std::cerr<<"DoenFr\n";
                        }
                    }
                }
            }
        }
    }

};

#endif
