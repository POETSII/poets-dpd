#ifndef naive_dpd_engine_half_merge_tbb_v3_hpp
#define naive_dpd_engine_half_merge_tbb_v3_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_merge_tbb.hpp"

class NaiveDPDEngineHalfMergeTBBV3
    : public NaiveDPDEngineHalfMergeTBB
{
public:
    std::vector<std::pair<double,double>> m_interaction_matrix;

    void Attach(WorldState *s) override
    {
        NaiveDPDEngineHalfMergeTBB::Attach(s);

        if(s){
            m_interaction_matrix.resize(s->interactions.size());
            for(unsigned i=0; i<m_interaction_matrix.size(); i++){
                m_interaction_matrix[i]={
                    s->interactions[i].conservative,
                    sqrt(s->interactions[i].dissipative)
                };
            }
        }
    }

private:
    virtual void step()
    {
        if(ForceLogging::logger()){
            step_impl<true>();
        }else{
            step_impl<false>();
        }
    }

    template<bool EnableLogging>
    void step_impl()
    {
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        if(EnableLogging && ForceLogging::logger()){
            ForceLogging::logger()->SetTime(m_state->t);
            ForceLogging::logger()->LogProperty("dt", 1, &m_state->dt);
            double seed_low=m_state->seed &0xFFFFFFFFul;
            double seed_high=m_state->seed>>32;
            ForceLogging::logger()->LogProperty("seed_lo", 1, &seed_low);
            ForceLogging::logger()->LogProperty("seed_high", 1, &seed_high);
            double t_hash_low=m_t_hash&0xFFFFFFFFul;
            double t_hash_high=m_t_hash>>32;
            ForceLogging::logger()->LogProperty("t_hash_lo", 1, &t_hash_low);
            ForceLogging::logger()->LogProperty("t_hash_high", 1, &t_hash_high);
            for(auto &b : m_state->beads){
                double h=b.get_hash_code().hash;
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(), "b_hash", 1, &h);
                double x[3]={b.x[0],b.x[1],b.x[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x",3,x);
                Bead bt=b;
                dpd_maths_core_half_step::update_mom(m_state->dt/2, bt);
                double v[3]={bt.v[0],bt.v[1],bt.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f",3,f);
            }
        }

        double dt=m_state->dt;

        // Move the beads, and then assign to cells based on x(t+dt)
        parallel_for_each_cell_blocked([&](Cell *pc){
            Cell &c=*pc;
            for(int bi=c.beads.size()-1; bi>=0; bi--){
                Bead *b=c.beads[bi];
                //std::cerr<<"In "<<c.pos<<" at "<<m_state->t<<"\n";
                dpd_maths_core_half_step::update_pos(dt, m_lengths, *b);
                unsigned index=world_pos_to_cell_index(b->x);
                if(index!=c.index){
                    //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                    auto &dst_cell=m_cells[index];
                    {
                        dst_cell.incoming_beads.push_back(b);
                    }
                    c.beads[bi]=c.beads.back();
                    c.beads.pop_back();
                }
            }
        });

        // At this point each cell will have most beads in c.beads, and might have some in c.incoming

        // Turns out more efficient to do explicitly than all at once
        parallel_for_each(m_cells, 1024, [&](Cell &c){
            transfer_incoming(c);
        });

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        parallel_for_each_cell_blocked([&](Cell *c){ process_cell<EnableLogging>(c); } );

        // Update all bonds
        parallel_for_each(m_non_monomers, m_non_monomer_grain, [&](const Polymer *p){
            const auto &pt = m_state->polymer_types.at(p->polymer_type);
            for(const auto &bond : pt.bonds){
                update_bond(*p, pt, bond);
            }
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(*p, pt, bond_pair);
            }
        });

        // Final mom
        parallel_for_each(m_state->beads, 1024, [&](Bead &b){
            dpd_maths_core_half_step::update_mom(m_state->dt, b);
        });

        if(EnableLogging && ForceLogging::logger()){
            for(auto &b : m_state->beads){

                if(ForceLogging::logger()){
                    double f[3]={b.f[0],b.f[1],b.f[2]};
                    ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_next",3,f);
                }
            }
        }

        m_state->t += 1;
    }

    // Idempotent function to moving incoming to beads. Should be fast in case where incoming is empty
    void transfer_incoming(Cell &c)
    {
        if(!c.incoming_beads.empty()){
            c.beads.insert(c.beads.end(), c.incoming_beads.begin(), c.incoming_beads.end());
            c.incoming_beads.clear();
        }
    };

    template<bool EnableLogging>
    void process_cell(Cell *c)
    {
        if(c->is_edge){
            process_cell_neighbours<EnableLogging,true>(*c);
        }else{
            process_cell_neighbours<EnableLogging,false>(*c);
        }
        update_intra_forces<EnableLogging>(*c);
    }

    template<bool EnableLogging,bool IsEdge>
    void process_cell_neighbours(Cell &c)
    {
        assert(c.neighbours.size()>=3);

        int i=c.neighbours.size()-1;
        Cell *next_next_next=c.neighbours[i-2];
        Cell *next_next=c.neighbours[i-1];
        Cell *next=c.neighbours[i];
        while(i>=0){
            Cell *curr=next;
            next=next_next;
            next_next=next_next_next;
            //Cell *curr=c.neighbours[i];
            --i;
            if(i>=2){
                next_next_next=c.neighbours[i-2];
                #ifdef __x86_64
                _mm_prefetch(&next_next_next->beads, _MM_HINT_T0);
                #endif
            }
            if(i>=1){
                #ifdef __x86_64
                _mm_prefetch(&next_next->beads[0], _MM_HINT_T0);
                #endif
            }
            if(next){
                #ifdef __x86_64
                for(auto p : next->beads){
                    _mm_prefetch(&p->x, _MM_HINT_T0);
                }
                #endif
            }
            update_inter_forces<EnableLogging,IsEdge>(c, *curr);
        }
    }

    template<bool EnableLogging>
    void update_bead_pair(Bead *hb, Bead *ob, const vec3r_t &dx, double dr2, vec3r_t &h_f)
    {
        assert(hb!=ob);

        double dr=pow_half(dr2);

        vec3r_t f;

        auto ia=m_interaction_matrix[ hb->bead_type*m_numBeadTypes+ob->bead_type ];
        
        dpd_maths_core_half_step::calc_force<EnableLogging,double,vec3r_t>(
            m_scaled_inv_root_dt,
            [&](unsigned a, unsigned b){ return ia.first; },
            [&](unsigned a, unsigned b){ return ia.second; },
            m_t_hash,
            dx, dr,
            0, 0, // kappa and r0
            *hb,
            *ob,
            f
        );

        h_f += f;
        ob->f -= f;
    }

};

#endif
