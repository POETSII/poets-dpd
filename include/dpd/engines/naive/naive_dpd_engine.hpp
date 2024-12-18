#ifndef naive_dpd_engine_hpp
#define naive_dpd_engine_hpp

#include "dpd/core/dpd_engine.hpp"

#include "dpd/maths/dpd_maths_core.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/dpd_state_validator.hpp"

#include <cassert>
#include <cmath>
#include <array>
#include <cfloat>
#include <iostream>

#include "dpd/core/logging.hpp"


/*
    The process to get from t to t+dt is:
    - UpdatePos : calculate position at time t+dt using values only from time t
        // Pre: x(t),v(t),f(t) all live
        x(t+dt) = x(t) + v(t)*dt + f(t)*(dt^2/2)
        v'(t+dt) = v(t) + f(t)*(lambda*dt),  where usually lambda=0.5
        // Post: x(t+dt),f(t),v(t),v'(t+dt) all live,   x(t) is dead
        // Note: v'(t+dt) is derived from v(t),f(t) in one mult_acc per dimension.

    - UpdateForce : calculate forces using positions at time t+dt
        // Pre: x(t+dt),f(t),v(t),v'(t+dt) all live
        f(t+dt) = (function of x(t+dt),v'(t+dt))
        // Post: f(t+dt),x(t+dt),f(t),v(t) all live

    - UpdateMom : calculate velocity at time t+dt
        v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2
        // Post: f(t+dt),v(t+dt),x(t+dt) all live, v(t) dead

    To keep state clean, the approach taken here is not to pre-calculate
    v'(t), and just to recalculate on the fly. This engine is already
    pretty inefficient anyway.
*/
template<bool UseMathsCore=false>
class NaiveDPDEngine
    : public DPDEngine
{
public:
    virtual double GetMaxBondLength() const
    {
        return 1000;
    }

    virtual void Attach(WorldState *state)
    {
        m_state=state;
        if(state){
            m_forces.resize(state->beads.size());
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        check_constraints_and_setup();

        for(unsigned i=0; i<nSteps; i++){
            step();
        }
    }

    bool CanSupportStationaryBeadTypes() const override
    { return true; }

    bool CanSupportSpatialDPDParameters() const override
    { return true; }

    bool CanSupportHookeanBonds() const override
    { return true; }

    bool CanSupportAngleBonds() const override
    { return true; }

    struct bp_hash
    {
        size_t operator()(const std::pair<const Bead *, const Bead *> ab) const
        {
            return std::hash<const Bead*>()(ab.first) + 19937 * std::hash<const Bead*>()(ab.second);
        }
    };
private:
    WorldState *m_state;

    struct Cell
    {
        unsigned index;
        vec3i_t pos;
        std::vector<Bead*> beads;

        // Spatially varying DPD only
        std::vector<InteractionStrength> local_interactions;
    };

    unsigned m_numBeadTypes;
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    vec3i_t m_dims; // Same as lengths, but integer
    std::vector<vec3r_t> m_forces;
    std::vector<Cell> m_cells;

    double m_inv_root_dt;

    uint64_t m_t_hash;

    bool m_has_spatially_varying_dpd_parameters=false;
    bool m_has_stationary_bead_types=false;
    

    std::unordered_set<std::pair<const Bead*,const Bead*>, bp_hash> m_seen_pairs;


    void check_constraints_and_setup()
    {
        validate(*m_state, 10);

        auto require=[](bool cond, const char *msg)
        {
            if(!cond){
                throw std::string(msg);
            }
        };

        require( m_state, "No state" );
        for(unsigned i=0; i<3; i++){
            require( round(m_state->box[i]) == m_state->box[i], "box must be integer aligned");
            require( m_state->box[i] >= 2, "Distance in each direction must be at least 2");
            m_lengths[i]=m_state->box[i];
            m_dims[i]=(int)m_state->box[i];
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_inv_root_dt=recip_pow_half(m_state->dt);

        m_has_stationary_bead_types=false;
        for(auto bt : m_state->bead_types){
            m_has_stationary_bead_types |= bt.stationary;
        }

        m_has_spatially_varying_dpd_parameters=false;
        for(auto i : m_state->interactions){
            if(!i.conservative.is_constant() || !i.dissipative.is_constant()){
                m_has_spatially_varying_dpd_parameters=true;
            }
        }

        std::unordered_map<std::string,double> bindings;

        // Reset all cell information
        m_cells.resize(calc_num_cells());
        unsigned ci=0;
        for(auto &c : m_cells){
            c.index=ci;
            c.pos=index_to_cell_pos(ci);
            assert(c.index==cell_pos_to_index(c.pos));
            c.beads.clear();

            if(!m_has_spatially_varying_dpd_parameters){
                c.local_interactions.clear();
            }else{
                c.local_interactions.resize(m_state->interactions.size());

                bindings["x"]=c.pos[0];
                bindings["y"]=c.pos[1];
                bindings["z"]=c.pos[2];
                bindings["ux"]= (c.pos[0] + 0.5 ) / m_dims[0] ;
                bindings["uy"]= (c.pos[1] + 0.5 ) / m_dims[1] ;
                bindings["uz"]= (c.pos[2] + 0.5 ) / m_dims[2] ;
                
                //std::cerr<<"  x="<<c.pos[0]<<" ";
                for(unsigned i=0; i<m_state->interactions.size(); i++){
                    c.local_interactions[i].conservative =
                        m_state->interactions[i].conservative.evaluate(bindings);
                    c.local_interactions[i].dissipative =
                        m_state->interactions[i].dissipative.evaluate(bindings);
                    //std::cerr<<" "<<c.local_interactions[i].conservative;
                    
                }
                //std::cerr<<"\n";
            }

            ++ci;
        }
    }

    unsigned calc_num_cells() const
    {
        unsigned n=1;
        for(unsigned i=0;i<3;i++){
            n *= m_dims[i];
        }
        return n;
    }

    vec3i_t world_pos_to_cell_pos(const vec3r_t &pos) const
    {
        vec3i_t res;
        for(unsigned i=0; i<3; i++){
            res[i] = (int)floor(pos[i]);
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
        return { (int)(index/(m_dims[1]*m_dims[2])) , (int)((index/m_dims[2])%m_dims[1]) , (int)(index%m_dims[2]) };
    }

    unsigned world_pos_to_cell_index(const vec3r_t &pos) const
    { return cell_pos_to_index(world_pos_to_cell_pos(pos)); }

    // This produces a pair (cell_loc,pos_adj)
    // - cell_loc : the cell position of the cell in direction dir
    // - delta : the delta to add to beads from the cell in direction dir to get correct relative position
    std::pair<vec3i_t,vec3r_t> make_relative_cell_pos(const vec3i_t &base, const vec3i_t &dir) const
    {
        vec3i_t pos;
        vec3r_t delta;
        for(unsigned i=0; i<3; i++){
            assert(-1 <= dir[i] && dir[i] <= +1);
            int raw = base[i] + dir[i];
            int base_is_left = raw < 0;
            int base_is_right = raw == m_lengths[i];

            pos[i] = raw + (base_is_left - base_is_right) * m_dims[i];

            delta[i] = (base_is_right - base_is_left) * m_lengths[i];
        }
        return {pos,delta};
    }

    vec3r_t calc_distance_from_to(const vec3r_t &base, const vec3r_t &other)
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

    void step()
    {
        m_t_hash = get_t_hash(m_state->t, m_state->seed);

        m_seen_pairs.clear();

        if(ForceLogging::logger()){
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
                double v[3]={b.v[0],b.v[1],b.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f",3,f);
            }
        }

        

        // Clear all cell information
        for(auto &c : m_cells){
            c.beads.clear();
        }

        m_forces.assign(m_state->beads.size(), vec3r_t{0,0,0});

        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &b : m_state->beads){
            update_bead_pos(&b);
            auto bpos=world_pos_to_cell_pos(b.x);
            auto bindex=cell_pos_to_index(bpos);
            auto &cell = m_cells.at( bindex );
            assert(cell.pos==bpos);
            cell.beads.push_back(&b);
            if(cell.beads.size() >= 16){
                std::cerr<<"  b="<<b.x<<", cell size="<<cell.beads.size()<<"\n";
            }

            if(ForceLogging::logger()){
                double x[3]={b.x[0],b.x[1],b.x[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x_next",3,x);
            }
        }

        // Calculate all the non-bonded forces
        // Each cell's force is calculated independently and over-written
        vec3i_t pos;
        for(pos[0]=0; pos[0]<m_lengths[0]; pos[0]++){
            for(pos[1]=0; pos[1]<m_lengths[1]; pos[1]++){
                for(pos[2]=0; pos[2]<m_lengths[2]; pos[2]++){
                    step_cell(pos);
                }
            }
        }

        for(auto &p : m_state->polymers){
            update_polymer_bond_forces(p);
        }

        for(auto &b : m_state->beads){
            double ff[3]={m_forces[b.bead_id][0], m_forces[b.bead_id][1], m_forces[b.bead_id][2]};
            if(ForceLogging::logger()){
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_now",3,ff);
            }
            update_bead_mom(&b);
            if(ForceLogging::logger()){
                double v[3]={b.v[0],b.v[1],b.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v_next",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f_next",3,f);
            }
        }

        m_state->t += 1;

        if(ForceLogging::logger()){
            ForceLogging::logger()->Flush();
        }
    }

    void step_cell(vec3i_t pos)
    {
        auto &home=m_cells[cell_pos_to_index(pos)];
        //std::cerr<<"  cell["<<pos<<"] -> "<<home.size()<<"\n";
        const InteractionStrength *home_interactions=nullptr;
        if(m_has_spatially_varying_dpd_parameters){
            home_interactions=&home.local_interactions[0];
        }

        vec3i_t dir;
        for(dir[0]=-1; dir[0]<=+1; dir[0]++){
            for(dir[1]=-1; dir[1]<=+1; dir[1]++){
                for(dir[2]=-1; dir[2]<=+1; dir[2]++){
                    auto [other_pos,other_delta] = make_relative_cell_pos(pos, dir);
                    const auto &other=m_cells[cell_pos_to_index(other_pos)];
                    assert(other.pos==other_pos);
                    const InteractionStrength *other_interactions=nullptr;
                    if(m_has_spatially_varying_dpd_parameters){
                        other_interactions=&other.local_interactions[0];
                    }

                    if(!home.beads.empty() && !other.beads.empty()){
                        //std::cerr<<"  base="<<pos<<", dir="<<dir<<", pos="<<other_pos<<", delta="<<other_delta<<"\n";
                        update_cell_forces(home.beads, home_interactions, other.beads, other_delta, other_interactions);
                    }
                }
            }
        }
    }

    void __attribute__((noinline)) update_cell_forces(
        std::vector<Bead*> &home, 
        const InteractionStrength *home_interactions,
        const std::vector<Bead*> &other,
        const vec3r_t &other_delta,
        const InteractionStrength *other_interactions
    )
    {
        for(Bead *hb : home){
            for(const Bead *ob : other)
            {
                update_bead_forces(hb, home_interactions, ob, other_delta, other_interactions);
            }
        }
    }


    static int32_t uint32_to_int32(uint32_t x)
    {
        // Sigh, Avoid undefined behaviour. Compiler should optimise it out.
        // https://stackoverflow.com/a/13208789
        int32_t res;
        if (x <= INT32_MAX) {
            res=static_cast<uint32_t>(x);
        }else{
            assert(x >= (uint32_t)INT32_MIN);
            res= static_cast<uint32_t>(x - INT32_MIN) + INT32_MIN;
        }
        assert(x == uint32_t(res) );  // int32_t -> uint32_t is well-defined
        return res;
    }

    double RandSym(const  BeadHash &s1, const BeadHash &s2) const
    { return RandSym(s1.hash, s2.hash); }

    double RandSym(uint32_t s1, uint32_t s2) const
    {
        uint32_t ru=hash_rng_sym(m_t_hash, s1,s2);
        int32_t rs=uint32_to_int32(ru);  // in [-2^31,2^31)
        //const double scale=ldexp(2.0,-32) / sqrt(1/3.0); // gives stddev of 1 (same as groot-warren paper)
        const double scale=ldexp(1.0, -32); // Gives range of [-0.5,0.5]  (same as Osprey-DPD)
        double u = rs * scale; 
        
        #if 0
        static double u_sum_sqr=0, u_sum=0;
        static unsigned u_count=0;
        u_sum_sqr += u*u;
        u_sum += u;
        u_count += 1;
        if((u_count % 10000)==0){
            std::cerr<<"ucount="<<u_count<<", mean="<<u_sum/u_count<<", std="<<sqrt(u_sum_sqr/u_count)<<"\n";
        }
        #endif
        
        //fprintf(stderr, "Naive : %llu, %u, %u -> %u, %f\n", m_t_hash, s1, s2, ru, u);

        return u;
    }

    void update_bead_forces(
        Bead *hb,
        const InteractionStrength *home_interactions,
        const Bead *ob,
        const vec3r_t &other_delta,
        const InteractionStrength *other_interactions
    )
    {
        if(hb==ob){
            return;
        }

        vec3r_t dx = vec3r_t(hb->x) - vec3r_t(ob->x) - other_delta; 
        double dr=dx.l2_norm();

        auto get_interaction_conservative=[&](unsigned atype, unsigned btype)
        {
            unsigned index=atype*m_numBeadTypes+btype;
            if(m_has_spatially_varying_dpd_parameters){
                double res=(home_interactions[index].conservative+other_interactions[index].conservative)*0.5;
                return res;
            }else{
                return (double)m_state->interactions[index].conservative;
            }
        };

        auto get_interaction_dissipative=[&](unsigned atype, unsigned btype)
        {
            unsigned index=atype*m_numBeadTypes+btype;
            if(m_has_spatially_varying_dpd_parameters){
                return (home_interactions[index].dissipative+other_interactions[index].dissipative)*0.5;
            }else{
                return (double)m_state->interactions[index].dissipative;
            }
        };


        if(UseMathsCore){
            vec3r_t f;
            if(dr < 1 && dr>=0.000000001){
                dpd_maths_core::calc_force(
                    (m_state->lambda * m_state->dt), pow_half( dpd_maths_core::kT * 24 / m_state->dt),
                    get_interaction_conservative,
                    get_interaction_dissipative,
                    m_t_hash,
                    dx, dr,
                    *hb,
                    *ob,
                    f
                );
                m_forces.at(hb->bead_id) += f;
            }
            return;
        }

        double inv_dr=1.0/dr;

        //std::cerr<<"ref : a="<<hb->get_bead_id()<<", b="<<ob->get_bead_id()<<", dx="<<dx<<"\n";


        auto rdx=calc_distance_from_to(hb->x, ob->x);
        double rdxr=rdx.l2_norm();
        //assert(rdxr <= 2*sqrt(2));

        if(dr>=1 || dr < 0.000000001){
            return;
        }

        double lambda_dt = (m_state->lambda * m_state->dt);
        vec3r_t hb_v = hb->v + hb->f * lambda_dt;
        vec3r_t ob_v = ob->v + ob->f * lambda_dt;
        
        vec3r_t dv = hb_v - ob_v;

        double wr = (1.0 - dr);
        double wr2 = wr*wr;
        
        double conForce = get_interaction_conservative(hb->bead_type,ob->bead_type)*wr;
        
        double rdotv = (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
		double gammap = get_interaction_dissipative(hb->bead_type,ob->bead_type)*wr2;

        double dissForce = -gammap*rdotv;
        double u = RandSym(hb->get_hash_code(), ob->get_hash_code());
        //fprintf(stderr, " r(%llu,%u,%u) -> %f\n", m_t_hash, hb->get_hash_code(), ob->get_hash_code(), u);
        double scaled_inv_root_dt=pow_half( dpd_maths_core::kT * 24 / m_state->dt);
        double randScale=pow_half(gammap) * scaled_inv_root_dt;
		double randForce = randScale * u;

        double scaled_force = (conForce + dissForce + randForce) * inv_dr;

        //std::cerr<<"  "<<hb->get_hash_code()<<" -> "<<ob->get_hash_code()<<" : "<<conForce<<", "<<dissForce<<", "<<randForce<<"\n";


        //std::cerr<<"  REF: home="<<hb->get_hash_code()<<", other="<<ob->get_hash_code()<<", t_hash="<<m_t_hash<<", dx="<<dx<<", r="<<dr<<", u="<<u<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<", hook=?\n";
        //std::cerr<<"     sqrt_gammap="<<sqrt(gammap)<<", rdotv="<<rdotv<<", sqrt(dissStrength)="<<sqrt(interactions.dissipative)<<"\n";

        vec3r_t f=dx * scaled_force;
        //std::cerr<<"ref :   dr="<<dr<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<"\n";
        m_forces.at(hb->bead_id) += f;
        /*if(hb->bead_id+1==2){
            std::cerr<<"  hb="<<hb->bead_id+1<<", ob="<<ob->bead_id+1<<", f_acc="<<m_forces.at(hb->bead_id)<<", f="<<f<<"\n";
        }*/

        if(ForceLogging::logger()){
            double ddx[3]={dx[0],dx[1],dx[2]};
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dx", 3,ddx);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dr", 1,&dr);
            double ddv[3]={dv[0],dv[1],dv[2]};
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dv", 3,ddv);
            double dd=get_interaction_dissipative(hb->bead_type,ob->bead_type);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-diss-strength", 1, &dd);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-invrootdt", 1, &scaled_inv_root_dt);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-gammap", 1, &gammap);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-rng", 1, &u);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-con", 1, &conForce);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-diss", 1,&dissForce);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-rng-scale",1, &randScale);
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"dpd-rand",1, &randForce);
            double ff[3]={f[0],f[1],f[2]};
            ForceLogging::logger()->LogBeadPairProperty(hb->get_hash_code(),ob->get_hash_code(),"f_next_dpd", 3,ff);
            
            if(!m_seen_pairs.insert({hb,ob}).second){
                fprintf(stderr, "  Seen %u %u twice\n", hb->get_hash_code().hash, ob->get_hash_code().hash);
            }
        }
    
        //std::cerr<<"Ref: t_hash="<<m_t_hash<<", h="<<hb->polymer_id<<", dx="<<rdx<<", dr="<<rdxr<<", f="<<f<<"\n";
    }

    void update_polymer_bond_forces(const Polymer &p)
    {
        const PolymerType &pt = m_state->polymer_types.at(p.polymer_type);
        for(const auto &bond : pt.bonds){
            update_hookean_bond(p, pt, bond);
        }

        for(const auto bi : p.bead_ids){
            const Bead *b = &m_state->beads[bi];
            //std::cerr<<"  bpo="<<b->get_hash_code()<<", fdpd="<<m_forces.at(b->bead_id)<<"\n";
        }

        for(const auto &bond_pair : pt.bond_pairs){
            update_angle_bond(p, pt, bond_pair);
        }
    }

    void update_hookean_bond(const Polymer &p, const PolymerType &, const Bond &b)
    {
        unsigned head_bead_index=p.bead_ids[b.bead_offset_head];
        unsigned tail_bead_index=p.bead_ids[b.bead_offset_tail];
        const Bead &head=m_state->beads.at(head_bead_index);
        const Bead &tail=m_state->beads.at(tail_bead_index);

        vec3r_t dx=calc_distance_from_to(head.x, tail.x);
        double r=dx.l2_norm();

        vec3r_t f;

        if(UseMathsCore){
            dpd_maths_core::calc_hookean_force<double,vec3r_t,vec3r_t>(
                b.kappa, b.r0, 
                dx, r, 1.0/r,
                f
            );
        }else{
            // The force scale is just kappa*(r-r0).
            double force=b.kappa*(r-b.r0);

            // Division by r is just to get (dx/r)
            f=dx * (force/r);
            //std::cerr<<" bond: r="<<r<<", force="<<force<<", fbond="<<f<<", dx="<<dx*(1/r)<<",  xh="<<head.x<<", xt="<<tail.x<<"\n";
        }
        
        m_forces[head.bead_id] += f;
        m_forces[tail.bead_id] -= f;
        
        if(ForceLogging::logger()){
            
            double ff[3]={f[0],f[1],f[2]};
            ForceLogging::logger()->LogBeadPairProperty(head.get_hash_code(), tail.get_hash_code(), "f_next_hookean", 3, ff);
            ForceLogging::logger()->LogBeadPairProperty(head.get_hash_code(), tail.get_hash_code(), "f_next_hookean_dx", 3, &dx[0]);
            ForceLogging::logger()->LogBeadPairProperty(head.get_hash_code(), tail.get_hash_code(), "f_next_hookean_dr", 1, &r);
            for(int i=0;i<3; i++){
                ff[i]=-ff[i];
            }
            double idx[3]={-dx[0],-dx[1],-dx[2]};
            ForceLogging::logger()->LogBeadPairProperty(tail.get_hash_code(), head.get_hash_code(), "f_next_hookean", 3, ff);
            ForceLogging::logger()->LogBeadPairProperty(tail.get_hash_code(), head.get_hash_code(), "f_next_hookean_dx", 3, idx);
            ForceLogging::logger()->LogBeadPairProperty(tail.get_hash_code(), head.get_hash_code(), "f_next_hookean_dr", 1, &r);
        }
    }

    void update_angle_bond(const Polymer &p, const PolymerType &pt, const BondPair &bp)
    {
        unsigned head_bond_index=bp.bond_offset_head;
        unsigned tail_bond_index=bp.bond_offset_tail;
        const Bond &head_bond=pt.bonds.at(head_bond_index);
        const Bond &tail_bond=pt.bonds.at(tail_bond_index);

        unsigned head_bead_index=head_bond.bead_offset_head;
        unsigned middle_bead_index=head_bond.bead_offset_tail;
        assert(middle_bead_index==tail_bond.bead_offset_head);
        unsigned tail_bead_index=tail_bond.bead_offset_tail;

        const Bead &head_bead=m_state->beads.at(p.bead_ids.at(head_bead_index));
        const Bead &middle_bead=m_state->beads.at(p.bead_ids.at(middle_bead_index));
        const Bead &tail_bead=m_state->beads.at(p.bead_ids.at(tail_bead_index));

        auto first=calc_distance_from_to(head_bead.x, middle_bead.x);
        auto second=calc_distance_from_to(middle_bead.x, tail_bead.x);

        //std::cerr<<m_state->t<<", ("<<bp.bond_offset_head<<","<<bp.bond_offset_tail<<"), first="<<first<<", second="<<second<<"\n";

        double FirstLength   = first.l2_norm();
        double SecondLength  = second.l2_norm();

        vec3r_t headForce, middleForce, tailForce;

        if(UseMathsCore){
            dpd_maths_core::calc_angle_force<double,vec3r_t,vec3r_t>(
                bp.kappa, cos(bp.theta0), sin(bp.theta0), 
                first, FirstLength,
                second, SecondLength,
                headForce, middleForce, tailForce
            );
        }else{
            const double magProduct = FirstLength*SecondLength;

            if(magProduct > 0.0001)
            {
                const double b1MagSq		= FirstLength*FirstLength;
                const double b2MagSq		= SecondLength*SecondLength;
                const double b1Dotb2		= first[0]*second[0] + first[1]*second[1] + first[2]*second[2];
                const double b1b2Overb1Sq	= b1Dotb2/b1MagSq;
                const double b1b2Overb2Sq	= b1Dotb2/b2MagSq;
                const double cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;
                const double Modulus=bp.kappa;

                double forceMag = 0.0;

                // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0
                if(absolute(b1Dotb2) > 0.000001)
                {
                    double Prefactor = sqrt(1.0/cosPhiSq - 1.0);
                    Prefactor = std::max(Prefactor, 0.000001);

                    // Add the restoring force depending on whether there is a preferred angle
                    // for the bond pair or not

                    double CosPhi0=cos(bp.theta0);

                    if(bp.theta0 > 0.0)
                    {
                        const double SinPhi0=sin(bp.theta0);
                        forceMag = Modulus*(CosPhi0 - SinPhi0/Prefactor)/magProduct;
                    }
                    else
                    {
                        forceMag = Modulus/magProduct;
                    }
                }
                else
                {
                    forceMag	= Modulus/magProduct;
                }

                //std::cerr<<"  Modulus="<<Modulus<<", phi0="<<bp.theta0<<", cosPhiSq="<<cosPhiSq<<", forceMag="<<forceMag<<", maxProduct="<<magProduct<<"\n";

                headForce =((first*b1b2Overb1Sq)-second) * forceMag;
                assert(isfinite(headForce));

                tailForce = (first- (second*b1b2Overb2Sq)) * forceMag;
                assert(isfinite(tailForce));

                middleForce = -(headForce + tailForce);
            }
        }

        //std::cerr<<"  dx01="<<first<<", dx12="<<second<<", r01="<<FirstLength<<", r12="<<SecondLength<<"\n";
        //std::cerr<<"Ref: fh="<<headForce<<", fm="<<middleForce<<", ft="<<tailForce<<"\n";

        //std::cerr<<"  bpo="<<head_bead.get_hash_code()<<", fdpd="<<m_forces[head_bead_index]<<", fangle="<<headForce<<", f="<<m_forces[head_bead_index]+headForce<<"\n";

        m_forces[head_bead.bead_id] += headForce;
        m_forces[tail_bead.bead_id] += tailForce;
        m_forces[middle_bead.bead_id] += middleForce;

        if(ForceLogging::logger()){
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_xhd", head_bead.x);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_xmd", middle_bead.x);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_xtl", tail_bead.x);
            
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_dx01", first);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_dx12", second);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_head", headForce);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_mid", middleForce);
            ForceLogging::logger()->LogBeadTripleProperty(head_bead.get_hash_code(), middle_bead.get_hash_code(), tail_bead.get_hash_code(), "f_next_angle_tail", tailForce);
        }

        if(m_state->t==5649){
            //std::cerr<<"Ref ; t="<<m_state->t<<", dx01="<<first<<", dx12="<<second<<"\n";
       }

        
    }

    // Pre: b.x==x(t), b.v==v(t), b.f==f(t)
    // Post: b.x==x(t+dt), b.v==v(t), b.f==f(t)
    void update_bead_pos(Bead *b) const
    {
        if(UseMathsCore){
            if(m_has_stationary_bead_types && m_state->bead_types[b->bead_type].stationary){
                // do nothing
            }else{
                dpd_maths_core::update_pos(
                    m_state->dt,
                    m_state->box,
                    *b
                );
            }
            return;
        }

        double dt=m_state->dt;

        //x(t+dt) = x(t) + v(t)*dt + a(t)/2*dt^2

        if(m_has_stationary_bead_types && m_state->bead_types[b->bead_type].stationary){
            // Do nothing
        }else{

            vec3r_t x = b->x + b->v*dt + b->f*(dt*dt/2);
            vec3r_t adj, xa;
            for(int i=0; i<3; i++){
                int add=x[i] < 0;
                int sub=x[i] >= m_state->box[i];
                adj[i] = (add-sub) * m_lengths[i];
                xa[i] = x[i] + adj[i];
                assert( 0 <= xa[i] && xa[i] <= m_state->box[i] );
                // Hack to get precise less than relationship on upper box.
                xa[i] = std::min( xa[i] , std::nexttoward(m_state->box[i], -DBL_MAX) );
                assert(xa[i] < m_state->box[i]);
            }

    #ifndef PDPD_TINSEL
            //std::cerr<<"  ref: x="<<xa<<", x'="<<b->x<<", v="<<b->v<<"\n";
    #endif

            b->x = xa;
            //std::cerr<<"Orig="<<orig<<", adj="<<adj<<", final="<<b->x<<"\n";
        }
    }

    void update_bead_mom(Bead *b)
    {
        if(UseMathsCore){
            if(m_has_stationary_bead_types && m_state->bead_types[b->bead_type].stationary){
                b->v.clear();
                b->f.clear();
            }else{
                dpd_maths_core::update_mom(
                    m_state->dt,
                    m_forces[b->bead_id],
                    *b
                );
            }
            return;
        }

        if(m_has_stationary_bead_types && m_state->bead_types[b->bead_type].stationary){
            b->v.clear();
            b->f.clear();
        }else{
            double dt=m_state->dt;

            // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

            vec3r_t fNow=m_forces[b->bead_id];
            m_forces[b->bead_id]=vec3r_t();

            b->v = b->v + (b->f + fNow) * (dt / 2);
            b->f = fNow;

            assert( isfinite(b->v) );
            assert( isfinite(b->f) );
        }
    }
};

#endif
