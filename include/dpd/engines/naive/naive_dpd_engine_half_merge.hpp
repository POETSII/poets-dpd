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
    static double now()
    {
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return ts.tv_sec+1e-9*ts.tv_nsec;
    };

public:
    friend class NaiveDPDEngineHalfMergeTBB;

    bool CanSupportStationaryBeadTypes() const override
    { return true; }

    bool CanSupportHookeanBonds() const override
    { return true; }

    bool CanSupportAngleBonds() const override
    { return true; }

    bool CanSupportSpatialDPDParameters() const override
    { return true; }

    std::string CanSupport(const WorldState *s) const override
    {
        for(int d=0;d<3;d++){
            int l=round(s->box[d]);
            require(l==s->box[d], "World bounds must be integer.");
            if( l%4 ){
                return "All dimensions must be a multiple of 4.";
            }
        }
        return DPDEngine::CanSupport(s);
    }

    virtual void Attach(WorldState *state)
    {
        m_state=state;
        if(m_state){
            m_timings=timings_t{};

            double start=now();
            check_constraints_and_setup();
            m_timings.configure += now()-start;
        }
    }

    virtual void sync_packed_to_beads()
    {
        for(Cell &c : m_cells){
            for(Packed &p : c.packed){
                auto b=p.bead;
                p.x=b->x;
                p.v=b->v;
                p.f=b->f;
            }
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        double start=now();

        for(unsigned i=0; i<nSteps; i++){
            step();
        }

        if(m_packed){
            sync_packed_to_beads();
        }

        double end=now();
        m_timings.execute_to_first_bead += end-start;
        m_timings.execute_to_last_bead += end-start;
    }

    virtual bool GetTimings(DPDEngine::timings_t &timings)
    {
        timings = m_timings;
        return true;
    }

protected:
    WorldState *m_state;
    timings_t m_timings;

    bool m_packed = false;

    struct interaction_strength
    {
        double conservative;
        double sqrt_dissipative;
    };

    struct Packed
    {
        uint32_t hash;
        vec3f_t x;
        vec3f_t v;
        vec3f_t f;
        Bead *bead; // Needed for angle bond lookups
    };

    struct Cell
    {
        unsigned index;
        vec3i_t pos;
        bool is_edge;
        std::vector<Cell*> neighbours;

        // Implementations use one or the other of raw beads or packed beads.
        std::vector<Bead*> beads;
        std::vector<Bead*> incoming_beads;

        std::vector<Packed> packed;
        std::vector<Packed> incoming_packed;

        std::vector<interaction_strength> local_interactions;
    };

    unsigned m_numBeadTypes;
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    vec3i_t m_dims;// Same as m_leNgths
    std::vector<Cell> m_cells;

    double m_scaled_inv_root_dt;
    uint64_t m_t_hash;

    std::vector<interaction_strength> m_global_interactions;

    bool m_has_spatial_dpd_interactions;
    bool m_has_stationary_bead_types;

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

        m_has_spatial_dpd_interactions=false;
        m_global_interactions.resize(m_state->interactions.size());
        for(unsigned index=0;index<m_state->interactions.size();index++){
            const auto &ii=m_state->interactions[index];
            if( !ii.conservative.is_constant() || !ii.dissipative.is_constant()){
                m_has_spatial_dpd_interactions=true;
                m_global_interactions.clear();
                break;
            }
            m_global_interactions[index].conservative=ii.conservative;
            m_global_interactions[index].sqrt_dissipative=sqrt(ii.dissipative);
        }

        m_has_stationary_bead_types=false;
        for(const auto &bt : m_state->bead_types){
            if(bt.stationary){
                std::cerr<<"Bead "<<bt.name<<" is stationary\n";
            }
            m_has_stationary_bead_types |= bt.stationary;
        }

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
            c.incoming_beads.clear();
            c.packed.clear();
            c.incoming_packed.clear();

            c.neighbours.clear();
            for(auto d : rel_nhood){
                vec3i_t neighbour=vec_wrap(c.pos+d, m_dims);
                c.neighbours.push_back(&m_cells.at(cell_pos_to_index(neighbour)));
            }

            if(m_has_spatial_dpd_interactions){
                std::unordered_map<std::string,double> bindings;
                bindings["x"]=c.pos[0];
                bindings["y"]=c.pos[1];
                bindings["z"]=c.pos[2];
                bindings["ux"]=(c.pos[0]+0.5) / m_state->box[0];
                bindings["uy"]=(c.pos[1]+0.5) / m_state->box[1];
                bindings["uz"]=(c.pos[2]+0.5) / m_state->box[2];
                c.local_interactions.resize(m_state->interactions.size());
                for(unsigned index=0; index<c.local_interactions.size(); index++){
                    c.local_interactions[index].conservative =
                        m_state->interactions[index].conservative.evaluate(bindings);
                    c.local_interactions[index].sqrt_dissipative =
                        sqrt(m_state->interactions[index].dissipative.evaluate(bindings));
                }
            }
        }

        if(!m_packed){
            for(auto &b : m_state->beads)
            {
                vec3i_t pos=floor(b.x);
                m_cells.at(cell_pos_to_index(pos)).beads.push_back(&b);
            }
        }else{
            for(auto &b : m_state->beads)
            {
                vec3i_t pos=floor(b.x);
                Packed p;
                p.bead=&b;
                p.f=b.f;
                p.v=b.v;
                p.x=b.x;
                p.hash=b.get_hash_code().hash;
                m_cells.at(cell_pos_to_index(pos)).packed.push_back(p);
            }
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
        assert(index<unsigned(m_dims[0]*m_dims[1]*m_dims[2]));
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
        assert(!m_packed);

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
                if(m_has_stationary_bead_types && m_state->bead_types[b->bead_type].stationary){
                    // do nothing
                    std::cerr<<"  not moving "<<b->polymer_id<<"\n";
                }else{
                    dpd_maths_core_half_step::update_pos(dt, m_lengths, *b);
                    unsigned index=world_pos_to_cell_index(b->x);
                    if(index!=c.index){
                        //std::cerr<<"Migrate at "<<m_state->t<<", "<<c.pos<<" -> "<<m_cells.at(index).pos<<"\n";
                        m_cells.at(index ).incoming_beads.push_back(b);
                        c.beads[bi]=c.beads.back();
                        c.beads.pop_back();
                    }
                }
            }
        }

        // At this point each cell will have most beads in c.beads, and might have some in c.incoming

        // Idempotent function to moving incoming to beads. Should be fast in case where incoming is empty
        auto transfer_incoming=[&](Cell &c)
        {
            if(!c.incoming_beads.empty()){
                c.beads.insert(c.beads.end(), c.incoming_beads.begin(), c.incoming_beads.end());
                c.incoming_beads.clear();
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
            if(m_has_stationary_bead_types && m_state->bead_types[b.bead_type].stationary){
                b.f.clear();
                b.v.clear();
            }else{
                dpd_maths_core_half_step::update_mom(m_state->dt, b);
            }
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
        assert(!m_packed);

        const interaction_strength *hinteractions = m_has_spatial_dpd_interactions
            ? &home.local_interactions[0]
            : &m_global_interactions[0];

        const interaction_strength *ointeractions = m_has_spatial_dpd_interactions
            ? &other.local_interactions[0]
            : nullptr;

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
            bool other_stationary = m_has_stationary_bead_types && m_state->bead_types[ob->bead_type].stationary;

            vec3r_t ob_x = vec3r_t(ob->x);
            if(IsEdge){
                ob_x += other_delta;
            }
            for(Bead *hb : home.beads){
                if(other_stationary && m_state->bead_types[hb->bead_type].stationary){
                    continue;
                }

                vec3r_t dx =  hb->x - ob_x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, hinteractions, ob, ointeractions, dx, dr2);
            }
        }
    }

    template<bool EnableLogging,bool IsEdge>
    void  __attribute__((noinline))  update_inter_forces_packed( Cell &home, Cell &other)
    {
        assert(m_packed);

        if(m_has_spatial_dpd_interactions){
            throw std::runtime_error("Not implemented for spatial interactions.");
        }

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

        for(Packed &ob : other.packed)
        {
            vec3r_t ob_x = vec3r_t(ob.x);
            if(IsEdge){
                ob_x += other_delta;
            }
            for(Packed &hb : home.packed){
                vec3r_t dx =  hb.x - ob_x;
                float dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<(float)MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, &m_global_interactions[0], ob, nullptr, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void  __attribute__((noinline))  update_intra_forces( Cell &home )
    {
        assert(!m_packed);

        const interaction_strength *interactions = m_has_spatial_dpd_interactions
            ? &home.local_interactions[0]
            : &m_global_interactions[0];

        for(int i=0; i<(int)home.beads.size()-1; i++){
            Bead *hb=home.beads[i];

            bool home_stationary = m_has_stationary_bead_types && m_state->bead_types[hb->bead_type].stationary;

            for(int j=i+1; j<(int)home.beads.size(); j++){
                Bead *ob=home.beads[j];

                if(home_stationary && m_state->bead_types[ob->bead_type].stationary){
                    continue;
                }

                vec3r_t dx =  hb->x - ob->x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, interactions, ob, interactions, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void  __attribute__((noinline))  update_intra_forces_packed( Cell &home )
    {
        assert(m_packed);
        if(m_has_spatial_dpd_interactions){
            throw std::runtime_error("Not implemented for spatial interactions.");
        }

        for(int i=0; i<(int)home.packed.size()-1; i++){
            Packed &hb=home.packed[i];
            for(int j=i+1; j<(int)home.packed.size(); j++){
                Packed &ob=home.packed[j];

                vec3f_t dx =  hb.x - ob.x;
                float dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<(float)MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, &m_global_interactions[0], ob, nullptr, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void update_bead_pair(Bead *hb, const interaction_strength *hinteractions, Bead *ob, const interaction_strength *ointeractions, const vec3r_t &dx, double dr2)
    {
        assert(hb!=ob);

        double dr=pow_half(dr2);

        vec3r_t f;

        auto get_conservative=[&](unsigned a, unsigned b)
        {
            unsigned index=a*m_numBeadTypes+b;
            return ointeractions
                ? (hinteractions[index].conservative+ointeractions[index].conservative)*0.5
                : hinteractions[index].conservative;
        };
        auto get_sqrt_dissipative=[&](unsigned a, unsigned b)
        {
            unsigned index=a*m_numBeadTypes+b;
            return ointeractions
                ? (hinteractions[index].sqrt_dissipative+ointeractions[index].sqrt_dissipative)*0.5
                : hinteractions[index].sqrt_dissipative;
        };
        
        dpd_maths_core_half_step::calc_force<EnableLogging,double,vec3r_t>(
            m_scaled_inv_root_dt,
            get_conservative,
            get_sqrt_dissipative,
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

    template<bool EnableLogging>
    void update_bead_pair_packed(Packed &hb, Packed &ob, const vec3f_t &dx, float dr2)
    {
        assert(hb.hash!=ob.hash);

        float dr=pow_half(dr2);

        vec3f_t f;
        
        dpd_maths_core_half_step::calc_force<EnableLogging,float,vec3f_t>(
            m_scaled_inv_root_dt,
            [&](unsigned a, unsigned b){ return m_state->interactions[a*m_numBeadTypes+b].conservative; },
            [&](unsigned a, unsigned b){ return pow_half(m_state->interactions[a*m_numBeadTypes+b].dissipative); },
            m_t_hash,
            dx, dr,
            0, 0, // kappa and r0
            hb,
            ob,
            f
        );

        hb.f += f;
        ob.f -= f;
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

        dpd_maths_core_half_step::calc_angle_force<double,vec3r_t,vec3r_t>(
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
