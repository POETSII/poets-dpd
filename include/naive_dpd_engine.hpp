#ifndef naive_dpd_engine_hpp
#define naive_dpd_engine_hpp

#include "dpd_engine.hpp"

#include "dpd_maths_core.hpp"

#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <array>


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

private:
    WorldState *m_state;

    unsigned m_numBeadTypes;
    vec3r_t m_origins; // origin in each dimension. Must be an integer, though encoded as real
    vec3r_t m_lengths; // dimensions in all directions Must be an integer, though encoded as real
    std::vector<vec3r_t> m_forces;
    std::vector<std::vector<Bead*>> m_cells;

    double m_inv_root_dt;

    uint64_t m_t_hash;


    void check_constraints_and_setup()
    {
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
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_inv_root_dt=1.0/sqrt(m_state->dt);
    }

    unsigned calc_num_cells() const
    {
        unsigned n=1;
        for(unsigned i=0;i<3;i++){
            n *= m_lengths[i];
        }
        return n;
    }

    vec3i_t world_pos_to_cell_pos(const vec3r_t &pos) const
    {
        vec3i_t res;
        for(unsigned i=0; i<3; i++){
            res[i] = floor(pos[i] - m_origins[i]);
            assert(0<=res[i] && res[i] < m_lengths[i]);
        }
        return res;
    }

    unsigned cell_pos_to_index(vec3i_t pos) const
    {
        unsigned index=0;
        unsigned last_dim=0;
        for(unsigned i=0; i<3; i++){
            assert(0<= pos[i] && pos[i] < m_lengths[i] );
            index = index * last_dim + pos[i];
            last_dim = m_lengths[i];
        }
        return index;
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

            pos[i] = raw + (base_is_left - base_is_right) * m_lengths[i];

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
        m_t_hash = next_t_hash(m_state->seed);

        // Clear all cell information
        m_cells.resize(calc_num_cells());
        for(auto &c : m_cells){
            c.clear();
        }

        // Move the beads, and then assign to cells based on x(t+dt)
        for(auto &b : m_state->beads){
            update_bead_pos(&b);
            m_cells.at( world_pos_to_cell_index(b.x) ).push_back(&b);
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
            update_bead_mom(&b);
        }

        m_state->t += m_state->dt;
    }

    void step_cell(vec3i_t pos)
    {
        auto &home=m_cells[cell_pos_to_index(pos)];
        //std::cerr<<"  cell["<<pos<<"] -> "<<home.size()<<"\n";

        vec3i_t dir;
        for(dir[0]=-1; dir[0]<=+1; dir[0]++){
            for(dir[1]=-1; dir[1]<=+1; dir[1]++){
                for(dir[2]=-1; dir[2]<=+1; dir[2]++){
                    auto [other_pos,other_delta] = make_relative_cell_pos(pos, dir);
                    const auto &other=m_cells[cell_pos_to_index(other_pos)];
                    if(!home.empty() && !other.empty()){
                        //std::cerr<<"  base="<<pos<<", dir="<<dir<<", pos="<<other_pos<<", delta="<<other_delta<<"\n";
                        update_cell_forces(home, other, other_delta);
                    }
                }
            }
        }
    }

    void update_cell_forces(std::vector<Bead*> &home, const std::vector<Bead*> &other, const vec3r_t &other_delta)
    {
        for(Bead *hb : home){
            for(const Bead *ob : other)
            {
                update_bead_forces(hb, ob, other_delta);
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

    double RandSym(uint32_t s1, uint32_t s2) const
    {
        uint32_t ru=hash_rng_sym(m_t_hash, s1,s2);
        int32_t rs=uint32_to_int32(ru);  // in [-2^31,2^31)
        //const double scale=ldexp(2.0,-32) / sqrt(1/3.0); // gives stddev of 1 (same as groot-warren paper)
        const double scale=ldexp(2.0, -32); // Gives range of [-0.5,0.5]  (same as Osprey-DPD)
        double u = rs * scale; 
        /*
        static double u_sum_sqr=0, u_sum=0;
        static unsigned u_count=0;
        u_sum_sqr += u*u;
        u_sum += u;
        u_count += 1;
        if((u_count % 10000)==0){
            std::cerr<<"ucount="<<u_count<<", mean="<<u_sum/u_count<<", std="<<sqrt(u_sum_sqr/u_count)<<"\n";
        }
        */
        return u;
    }

    void update_bead_forces(Bead *hb, const Bead *ob, const vec3r_t &other_delta)
    {
        if(hb==ob){
            return;
        }

        vec3r_t dx = vec3r_t(hb->x) - vec3r_t(ob->x) - other_delta; 
        double dr=dx.l2_norm();

        if(UseMathsCore){
            vec3r_t f;
            if(dr < 1 && dr>=0.000000001){
                dpd_maths_core::calc_force(
                    (m_state->lambda * m_state->dt), (1.0/sqrt(m_state->dt)),
                    [&](unsigned a, unsigned b){ return m_state->interactions[a*m_numBeadTypes+b].conservative; },
                    [&](unsigned a, unsigned b){ return m_state->interactions[a*m_numBeadTypes+b].dissipative; },
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
        
        const auto &interactions=m_state->interactions[hb->bead_type*m_numBeadTypes+ob->bead_type];

        double lambda_dt = (m_state->lambda * m_state->dt);
        vec3r_t hb_v = hb->v + hb->f * lambda_dt;
        vec3r_t ob_v = ob->v + ob->f * lambda_dt;
        
        vec3r_t dv = hb_v - ob_v;

        double wr = (1.0 - dr);
        double wr2 = wr*wr;
        
        double conForce = interactions.conservative*wr;
        
        double rdotv = (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
		double gammap = interactions.dissipative*wr2;

        double dissForce = -gammap*rdotv;
        double u = RandSym(hb->get_hash_code(), ob->get_hash_code());
		double randForce = sqrt(gammap) * (1.0/sqrt(m_state->dt)) * u;

        double scaled_force = (conForce + dissForce + randForce) * inv_dr;

        //std::cerr<<"  "<<hb->get_hash_code()<<" -> "<<ob->get_hash_code()<<" : "<<conForce<<", "<<dissForce<<", "<<randForce<<"\n";


        //std::cerr<<"  REF: home="<<hb->get_hash_code()<<", other="<<ob->get_hash_code()<<", t_hash="<<m_t_hash<<", dx="<<dx<<", r="<<dr<<", u="<<u<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<", hook=?\n";
        //std::cerr<<"     sqrt_gammap="<<sqrt(gammap)<<", rdotv="<<rdotv<<", sqrt(dissStrength)="<<sqrt(interactions.dissipative)<<"\n";

        vec3r_t f=dx * scaled_force;
        //std::cerr<<"ref :   dr="<<dr<<", con="<<conForce<<", diss="<<dissForce<<", ran="<<randForce<<"\n";
        m_forces.at(hb->bead_id) += f;
        //std::cerr<<"  r="<<dr<<", force="<<scaled_force*dr<<", f="<<m_forces.at(hb->bead_id)<<"\n";
    
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
            dpd_maths_core::calc_hookean_force(
                b.kappa, b.r0, 
                dx, r, 1.0/r,
                f
            );
        }else{
            // The force scale is just kappa*(r-r0).
            double force=b.kappa*(r-b.r0);

            // Division by r is just to get (dx/r)
            f=dx * (force/r);
            //std::cerr<<" bond: r="<<r<<", force="<<force<<", fbond="<<f<<"\n";
        }
        
        m_forces[head_bead_index] += f;
        m_forces[tail_bead_index] -= f;
    }

    void update_angle_bond(const Polymer &p, const PolymerType &pt, const BondPair &bp)
    {
        unsigned head_bond_index=p.bead_ids[bp.bond_offset_head];
        unsigned tail_bond_index=p.bead_ids[bp.bond_offset_tail];
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

        double FirstLength   = first.l2_norm();
        double SecondLength  = second.l2_norm();

        vec3r_t headForce, middleForce, tailForce;

        if(UseMathsCore){
            dpd_maths_core::calc_angle_force(
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
                if(fabs(b1Dotb2) > 0.000001)
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

        m_forces[head_bead_index] += headForce;
        m_forces[tail_bead_index] += tailForce;
        m_forces[middle_bead_index] += middleForce;

        
    }

    // Pre: b.x==x(t), b.v==v(t), b.f==f(t)
    // Post: b.x==x(t+dt), b.v==v(t), b.f==f(t)
    void update_bead_pos(Bead *b) const
    {
        if(UseMathsCore){
            dpd_maths_core::update_pos(
                m_state->dt,
                m_state->box,
                *b
            );
            return;
        }

        double dt=m_state->dt;

        //x(t+dt) = x(t) + v(t)*dt + a(t)/2*dt^2

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
        }
        b->x = xa;
        //std::cerr<<"Orig="<<orig<<", adj="<<adj<<", final="<<b->x<<"\n";
    }

    void update_bead_mom(Bead *b)
    {
        if(UseMathsCore){
            dpd_maths_core::update_mom(
                m_state->dt,
                m_forces[b->bead_id],
                *b
            );
            return;
        }

        double dt=m_state->dt;

       // v(t+dt) = v(t) + dt*(f(t)+f(t+dt))/2

       vec3r_t fNow=m_forces[b->bead_id];
       m_forces[b->bead_id]=vec3r_t();

       b->v = b->v + (b->f + fNow) * (dt / 2);
       b->f = fNow;

       assert( isfinite(b->v) );
       assert( isfinite(b->f) );
    }
};

#endif
