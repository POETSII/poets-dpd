#ifndef naive_dpd_engine_half_merge_avx_hpp
#define naive_dpd_engine_half_merge_avx_hpp

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
#include <array>

#include <immintrin.h>

class NaiveDPDEngineHalfMergeAVX
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

    std::string CanSupport(const WorldState *s) const override
    {
        for(const PolymerType &pt : s->polymer_types){
            if(pt.bead_types.size()>1){
                return "Polymersnot yet supported.";
            }
        }

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
            m_timings=timings_t{};

            double start=now();
            check_constraints_and_setup();
            m_timings.configure += now()-start;
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        double start=now();

        for(unsigned i=0; i<nSteps; i++){
            step();
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

    struct Packed
    {
        uint32_t hash;
        vec3f_t x;
        vec3f_t v;
        vec3f_t f;
    };

    const unsigned VEC_WIDTH = 4;

    #error "Another half-arsed attempt"

    struct Cell
    {
        unsigned index;
        vec3i_t pos;
        bool is_edge;
        std::array<uint32_t,13> neighbour_indices;

        uint32_t n; // Number of local cells
        uint32_t todo; // Number of local cells not yet moved. 0 <= todo <= n

        float vec_x[3][VEC_WIDTH];
        float vec_v[3][VEC_WIDTH];
        uint32_t vec_hash[VEC_WIDTH];
        float vec_f[3][VEC_WIDTH];

        std::vector<Packed> packed;

        void remove(unsigned i)
        {
            assert(i<=n);
            if(i>=VEC_WIDTH){
                // Bead to deleteis in packed area
                if( i<n ){
                   packed[i-VEC_WIDTH] = packed.back();
                }
                packed.pop_back();
            }else if(i==VEC_WIDTH-1){
                // Bead to delete is the last beack in the vec area
                if(i<n){
                    // ... and there are also beads in the packed area
                    for(int d=0; d<3; d++){
                        vec_x[d][i] = packed.back().x[d];
                        vec_v[d][i] = packed.back().v[d];
                        vec_f[d][i] = packed.back().f[d];
                    }
                    vec_hash[i] = packed.back().hash;
                    packed.pop_back();
                }
            }else{
                // Bead to delete is in vec area and not the last
                if(i<n){
                    // .. and there is another vec bead after
                    for(int d=0; d<3; d++){
                        vec_x[d][i] = vec_x[d][n-1];
                        vec_v[d][i] = vec_v[d][n-1];
                        vec_f[d][i] = vec_f[d][n-1];
                    }
                    vec_hash[i] = vec_hash[n-1];
                }
            }
            --n;
        }

        void insert(const Packed &b)
        {
            if(n>=VEC_WIDTH){
                packed.push_back(b);
            }else{
                for(int d=0; d<3; d++){
                    vec_x[d][n] = b.x[d];
                    vec_v[d][n] = b.v[d];
                    vec_f[d][n] = b.f[d];
                }
                vec_hash[d][n] = b.hash[d];
            }
            ++n;
        }
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

            c.n=0;
            c.todo=0;
            c.packed.clear();

            c.neighbours.clear();
            unsigned offset=0;
            for(auto d : rel_nhood){
                vec3i_t neighbour=vec_wrap(c.pos+d, m_dims);
                c.
                c.neighbours[offset++]=cell_pos_to_index(neighbour);
            }
        }

        for(auto &b : m_state->beads)
        {
            vec3i_t pos=floor(b.x);
            Packed p;
            p.f=b->f;
            p.v=b->v;
            p.x=b->x;
            p.hash=b->GetHash();
            m_cells.at(cell_pos_to_index(pos)).insert(p);
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
                    m_cells.at(index ).incoming_beads.push_back(b);
                    c.beads[bi]=c.beads.back();
                    c.beads.pop_back();
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
        assert(!m_packed);

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

    template<bool EnableLogging,bool IsEdge>
    void  __attribute__((noinline))  update_inter_forces_packed( Cell &home, Cell &other)
    {
        assert(m_packed);

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

                update_bead_pair<EnableLogging>(hb, ob, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void  __attribute__((noinline))  update_intra_forces( Cell &home )
    {
        assert(!m_packed);

        for(int i=0; i<(int)home.beads.size()-1; i++){
            Bead *hb=home.beads[i];
            for(int j=i+1; j<(int)home.beads.size(); j++){
                Bead *ob=home.beads[j];

                vec3f_t dx =  hb->x - ob->x;
                double dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<MIN_DISTANCE_CUTOFF_SQR){
                    continue;
                }

                update_bead_pair<EnableLogging>(hb, ob, dx, dr2);
            }
        }
    }

    template<bool EnableLogging>
    void  __attribute__((noinline))  update_intra_forces_packed( Cell &home )
    {
        assert(m_packed);

        for(int i=0; i<(int)home.packed.size()-1; i++){
            Packed &hb=home.packed[i];
            for(int j=i+1; j<(int)home.packed.size(); j++){
                Packed &ob=home.packed[j];

                vec3f_t dx =  hb->x - ob->x;
                float dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if(dr2 >= 1 || dr2<(float)MIN_DISTANCE_CUTOFF_SQR){
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
            *hb,
            *ob,
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
