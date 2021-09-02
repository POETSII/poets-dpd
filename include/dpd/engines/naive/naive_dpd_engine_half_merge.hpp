#ifndef naive_dpd_engine_half_merge_hpp
#define naive_dpd_engine_half_merge_hpp

#include "dpd/core/dpd_engine.hpp"
#include "dpd/engines/naive/naive_dpd_engine.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"
#include "dpd/core/make_nhood.hpp"

#include <cassert>
#include <cmath>
#include <array>
#include <unordered_set>
#include <mutex>

#include <immintrin.h>

class NaiveDPDEngineHalfMerge
    : public DPDEngine
{
public:
    friend class NaiveDPDEngineHalfMergeTBB;

    std::string CanSupport(const WorldState *s) const override
    {
        for(int d=0;d<3;d++){
            int l=round(s->box[d]);
            require(l==s->box[d], "World bounds must be integer.");
            if( l%4 ){
                return "All dimensions must be a multiple of 4.";
            }
        }
        return {};
    }

    virtual void Attach(WorldState *state)
    {
        m_state=state;
        if(m_state){
            check_constraints_and_setup();
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        for(unsigned i=0; i<nSteps; i++){
            step();
        }
    }

private:
    WorldState *m_state;

    struct Cell
    {
        unsigned index;
        vec3i_t pos;
        bool is_edge;
        std::vector<Bead*> beads;
        std::vector<Cell*> neighbours;
        std::vector<Bead*> incoming;

        std::mutex mutex;
    };

    unsigned m_numBeadTypes;
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    vec3i_t m_dims;// Same as m_leNgths
    std::vector<Cell> m_cells;

    double m_scaled_inv_root_dt;
    uint64_t m_t_hash;

    void require(bool cond, const char *msg) const
    {
        if(!cond){
            throw std::string(msg);
        }
    };

    void check_constraints_and_setup()
    {
        require( m_state, "No state" );
        for(unsigned i=0; i<3; i++){
            require( round(m_state->box[i]) == m_state->box[i], "box must be integer aligned");
            require( m_state->box[i] >= 2, "Distance in each direction must be at least 2");
            m_lengths[i]=m_state->box[i];
            m_dims[i]=(int)m_state->box[i];
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_scaled_inv_root_dt=pow_half(24*dpd_maths_core_half_step::kT / m_state->dt);

        std::vector<vec3i_t> rel_nhood=make_relative_nhood_forwards(/*excludeCentre*/true);

        m_cells=std::vector<Cell>( m_dims[0]*m_dims[1]*m_dims[2] );
        for(unsigned i=0; i<m_cells.size(); i++){
            auto &c = m_cells[i];
            c.index=i;
            c.pos=index_to_cell_pos(i);
            assert(c.index==cell_pos_to_index(c.pos));
            c.is_edge=false;
            for(int d=0; d<3; d++){
                c.is_edge |= c.pos[d]==0 || c.pos[d]==m_dims[d]-1;
            }
            c.beads.clear();
            c.incoming.clear();

            c.neighbours.clear();
            for(auto d : rel_nhood){
                vec3i_t neighbour=vec_wrap(c.pos+d, m_dims);
                c.neighbours.push_back(&m_cells.at(cell_pos_to_index(neighbour)));
            }
        }

        for(auto &b : m_state->beads)
        {
            vec3i_t pos=floor(b.x);
            m_cells.at(cell_pos_to_index(pos)).beads.push_back(&b);
        }
    }

    vec3i_t world_pos_to_cell_pos(const vec3r_t &pos) const
    {
        vec3i_t res=floor(pos);
        for(unsigned i=0; i<3; i++){
            assert(0<=res[i] && res[i] < m_dims[i]);
        }
        return res;
    }

    unsigned cell_pos_to_index(vec3i_t pos) const
    {
        return pos[0]*m_dims[1]*m_dims[2] + pos[1] * m_dims[2] + pos[2];
    }

    vec3i_t index_to_cell_pos(unsigned index) const
    {
        assert(index<m_dims[0]*m_dims[1]*m_dims[2]);
        return { (int)(index/(m_dims[1]*m_dims[2])) , int((index/m_dims[2])%m_dims[1]) , int(index%m_dims[2]) };
    }

    unsigned world_pos_to_cell_index(const vec3r_t &pos) const
    { return cell_pos_to_index(world_pos_to_cell_pos(pos)); }

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
        for(auto &c : m_cells){
            for(int bi=c.beads.size()-1; bi>=0; bi--){
                Bead *b=c.beads[bi];
                //std::cerr<<"In "<<c.pos<<" at "<<m_state->t<<"\n";
                dpd_maths_core_half_step::update_pos(dt, m_lengths, *b);
                unsigned index=world_pos_to_cell_index(b->x);
                if(index!=c.index){
                    //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                    m_cells.at(index ).incoming.push_back(b);
                    c.beads[bi]=c.beads.back();
                    c.beads.pop_back();
                }
            }
        }

        // At this point each cell will have most beads in c.beads, and might have some in c.incoming

        // Idempotent function to moving incoming to beads. Should be fast in case where incoming is empty
        auto transfer_incoming=[&](Cell &c)
        {
            if(!c.incoming.empty()){
                c.beads.insert(c.beads.end(), c.incoming.begin(), c.incoming.end());
                c.incoming.clear();
            }
        };

        // Calculate all the DPD and 2-bead bond forces
        // Each cell's force is calculated independently
        for(auto &c : m_cells){
            transfer_incoming(c);
            update_intra_forces<EnableLogging>(c);
            if(c.is_edge){
                for(auto *n : c.neighbours){
                    transfer_incoming(*n);
                    update_inter_forces<EnableLogging,true>(c, *n);
                }
            }else{
                for(auto *n : c.neighbours){
                    transfer_incoming(*n);
                    update_inter_forces<EnableLogging,false>(c, *n);
                }
            }
        }

        // Update all bonds
        for(const auto &p : m_state->polymers){
            const auto &pt = m_state->polymer_types.at(p.polymer_type);
            for(const auto &bond : pt.bonds){
                update_bond(p, pt, bond);
            }
            for(const auto &bond_pair : pt.bond_pairs){
                update_angle_bond(p, pt, bond_pair);
            }
        }

        // Final mom
        for(auto &b : m_state->beads){
            dpd_maths_core_half_step::update_mom(m_state->dt, b);
        }

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

    template<bool EnableLogging,bool IsEdge>
    void  __attribute__((noinline))  update_inter_forces( Cell &home, Cell &other)
    {
        vec3r_t other_delta;
        if(IsEdge){
            for(int d=0; d<3; d++){
                if(home.pos[d]==0){
                    if(other.pos[d]==m_dims[d]-1){
                        other_delta[d]=-m_dims[d];
                    }
                }else if(home.pos[d]==m_dims[d]-1){
                    if(other.pos[d]==0){
                        other_delta[d]=m_dims[d];
                    }
                }
            }
        }

        for(Bead *ob : other.beads)
        {
            vec3r_t ob_x = vec3r_t(ob->x);
            if(IsEdge){
                ob_x += other_delta;
            }
            for(Bead *hb : home.beads){
                vec3r_t dx =  hb->x - ob_x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, ob, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void  __attribute__((noinline))  update_intra_forces( Cell &home )
    {
        for(int i=0; i<(int)home.beads.size()-1; i++){
            Bead *hb=home.beads[i];
            for(int j=i+1; j<(int)home.beads.size(); j++){
                Bead *ob=home.beads[j];

                vec3r_t dx =  hb->x - ob->x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, ob, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void update_bead_pair(Bead *hb, Bead *ob, const vec3r_t &dx, double dr2)
    {
        assert(hb!=ob);

        double dr=pow_half(dr2);

        vec3r_t f;
        
        dpd_maths_core_half_step::calc_force<EnableLogging,double,vec3r_t>(
            m_scaled_inv_root_dt,
            [&](unsigned a, unsigned b){ return m_state->interactions[a*m_numBeadTypes+b].conservative; },
            [&](unsigned a, unsigned b){ return pow_half(m_state->interactions[a*m_numBeadTypes+b].dissipative); },
            m_t_hash,
            dx, dr,
            0, 0, // kappa and r0
            *hb,
            *ob,
            f
        );

        hb->f += f;
        ob->f -= f;
    }

    vec3r_t calc_distance_from_to(const vec3r_t &base, const vec3r_t &other) const
    {
        vec3r_t res;
        for(unsigned i=0; i<3; i++){
            double len=m_lengths[i];
            double dx=other[i]-base[i];
            int wrap_neg = dx > len*0.5;
            int wrap_pos = dx < -len*0.5;
            res[i] = dx + (wrap_pos - wrap_neg) * len;
        }
        return res;
    }

    void update_bond(const Polymer &p, const PolymerType &pt, const Bond &b) const
    {
        unsigned head_index=p.bead_ids.at(b.bead_offset_head);
        unsigned tail_index=p.bead_ids.at(b.bead_offset_tail);
        Bead &head_bead=m_state->beads[head_index];
        Bead &tail_bead=m_state->beads[tail_index];

        vec3r_t dx=calc_distance_from_to(head_bead.x, tail_bead.x);
        double dr=dx.l2_norm();

        double dr0=b.r0-dr;
        double hookeanForce=b.kappa*dr0;

        vec3r_t f=dx * (hookeanForce/dr);
        head_bead.f -= f;
        tail_bead.f += f;
    }

    void update_angle_bond(const Polymer &p, const PolymerType &pt, const BondPair &bp) const
    {
        unsigned head_bond_index=bp.bond_offset_head;
        unsigned tail_bond_index=bp.bond_offset_tail;
        const Bond &head_bond=pt.bonds.at(head_bond_index);
        const Bond &tail_bond=pt.bonds.at(tail_bond_index);

        unsigned head_bead_index=head_bond.bead_offset_head;
        unsigned middle_bead_index=head_bond.bead_offset_tail;
        assert(middle_bead_index==tail_bond.bead_offset_head);
        unsigned tail_bead_index=tail_bond.bead_offset_tail;

        Bead &head_bead=m_state->beads.at(p.bead_ids.at(head_bead_index));
        Bead &middle_bead=m_state->beads.at(p.bead_ids.at(middle_bead_index));
        Bead &tail_bead=m_state->beads.at(p.bead_ids.at(tail_bead_index));

        auto first=calc_distance_from_to(head_bead.x, middle_bead.x);
        auto second=calc_distance_from_to(middle_bead.x, tail_bead.x);

        double FirstLength   = first.l2_norm();
        double SecondLength  = second.l2_norm();

        vec3r_t headForce, middleForce, tailForce;

        dpd_maths_core_half_step::calc_angle_force(
            bp.kappa, cos(bp.theta0), sin(bp.theta0), 
            first, FirstLength,
            second, SecondLength,
            headForce, middleForce, tailForce
        );

        head_bead.f += headForce;
        tail_bead.f += tailForce;
        middle_bead.f += middleForce;
    }
};

#endif
