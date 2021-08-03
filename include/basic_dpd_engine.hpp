#ifndef basic_dpd_engine_hpp
#define basic_dpd_engine_hpp

#include "dpd_engine.hpp"

#include "dpd_maths_core_half_step.hpp"

#include "vec3.hpp"
#include "make_nhood.hpp"

#include <cassert>
#include <cmath>
#include <array>
#include <cstdint>



/*
    This is a minimal engine which models the data requirements of working
    in POETS. So it is split into cells, where each call contains all
    information needed for the beads in that cell. When beads are moved
    they also contain all information needed to manage bonds, so bond
    information is embedded directly in the beads.

    This means there are the following restrictions:
    - Each bead can have at most 3 bonded pairs
    - Each bead can have at most 2 angle bonds

    Based on looking at input files from Julian, there are also the following
    conveniences:
    - All bonds have the same kappa (128) and r0 (0.5)
    - All angle bonds are straight
    - There are at most two bead kappas
*/
class BasicDPDEngine
    : public DPDEngine
{
public:
    static constexpr int MAX_BONDS_PER_BEAD=3;
    static constexpr int MAX_ANGLE_BONDS_PER_BEAD=1;
private:

    friend class BasicDPDEngineV2;
    friend class BasicDPDEngineV3;
    friend class BasicDPDEngineV4;
    friend class BasicDPDEngineV3Raw;
    friend class BasicDPDEngineV4Raw;
    friend class BasicDPDEngineV5Raw;

    static uint32_t make_bead_id(bool is_monomer, unsigned polymer_id, unsigned polymer_offset, unsigned bead_type)
    {
        uint32_t bead_id;
        if(is_monomer){
            assert(polymer_offset==0);
            bead_id = 0x8000'0000 | (polymer_id<<4) | bead_type;
        }else{
            bead_id = (polymer_offset<<24) | (polymer_id<<4) | bead_type;
        }
        return bead_id;
    }
    struct bead_id_t
    {
        /*  1ppp pppp pppp pppp  pppp pppp pppp bbbb
            0ooo oooo pppp pppp  pppp pppp pppp bbbb
            p = polymer id
            o = polymer offset
            b = bead type

            A polymer id plus offset is just shifted right by 4
            0000 1ppp pppp pppp  pppp pppp pppp pppp
            0000 0ooo oooo pppp  pppp pppp pppp pppp
            This is generally cheaper than masking
        */
        uint32_t bead_id;

        void operator=(uint32_t x)
        {
            bead_id=x;
        }

        void set_bead_id(bool _is_monomer, uint32_t polymer_id, uint32_t polymer_offset, uint32_t bead_type)
        {
            bead_id=make_bead_id(_is_monomer, polymer_id, polymer_offset, bead_type);
            assert(is_monomer() == _is_monomer);
            assert(polymer_id == get_polymer_id());
            assert(is_monomer() ? polymer_offset==0 : polymer_offset == get_polymer_offset());
            assert(bead_type==get_bead_type());
        }

        uint32_t get_bead_type() const
        { return bead_id & 0xF; }

        bool is_monomer() const
        { return bead_id>>31; }

        uint32_t get_polymer_id() const
        { return (bead_id>>4) & (is_monomer() ? 0x07FF'FFFFul : 0x000F'FFFFul); }

        uint32_t get_polymer_offset() const
        { return is_monomer() ? 0 : (bead_id>>24)&0x7F; }

        uint32_t make_hash_from_offset(unsigned offset)
        {
            assert(!is_monomer());
            return ((bead_id>>4)&0x000F'FFFFul) | (offset<<20);
        }

        uint32_t get_hash_code() const
        {
            return bead_id>>4;
        }

    };
    static_assert(sizeof(bead_id_t)==4);

    struct bead_view_t
    {
        bead_id_t id;
        vec3f_t x;
        vec3f_t v;

        uint32_t get_hash_code() const
        { return id.get_hash_code(); }

        uint32_t get_bead_type() const
        { return id.get_bead_type(); }
    };
    static_assert(sizeof(bead_view_t)==28);

    struct angle_bond_info_t
    {
        uint8_t partner_head;
        uint8_t partner_tail;
        uint8_t kappa;   // Kappa is just an integer. Typically this is quite small, e.g. 5 or 15. Should be less than 255
        uint8_t _pad_;
    };
    static_assert(sizeof(angle_bond_info_t)==4);

    struct bead_exchange_t
        : bead_view_t
    {
        vec3f_t f;

        uint8_t bond_partners[MAX_BONDS_PER_BEAD+1]; // -1 means no bond
        static_assert(sizeof(bond_partners)==4); // Keep aligned

        angle_bond_info_t angle_bonds[MAX_ANGLE_BONDS_PER_BEAD];
    };
    static_assert(sizeof(bead_exchange_t)==48);

    struct bead_resident_t
        : bead_exchange_t
    {
        // Nothing here at the moment
    };
    static_assert(sizeof(bead_resident_t)==48);

    struct force_input_t
    {
        uint32_t target_hash;
        vec3f_t f;
    };

    void get_bond_info(
        const bead_resident_t &home,
        const bead_id_t &other_id,
        int &angle_partner_count,
        float &kappa,
        float &r0
    ){
        angle_partner_count=false;
        kappa=0.0f;
        r0=m_bond_r0;
        if(home.id.is_monomer()) return;
        if(other_id.is_monomer()) return;
        if(home.id.get_polymer_id() != other_id.get_polymer_id()) return;
        auto other_polymer_offset=other_id.get_polymer_offset();
        for(unsigned i=0; i<MAX_BONDS_PER_BEAD; i++){
            if(0xFF == home.bond_partners[i]){
                break;
            }
            if(other_polymer_offset==home.bond_partners[i]){
                kappa=m_bond_kappa;
                break;
            }
        }
        if(kappa==0.0f){
            return;
        }
        for(unsigned i=0; i<MAX_ANGLE_BONDS_PER_BEAD; i++){
            if(home.angle_bonds[i].partner_head==0xFF){
                return;
            }
            if(home.angle_bonds[i].partner_head==other_polymer_offset){
                angle_partner_count++;
            }
            if(home.angle_bonds[i].partner_tail==other_polymer_offset){
                angle_partner_count++;
            }
        }
    }

public:
    std::string CanSupport(const WorldState *state) const override
    {
        float bond_kappa=nanf("");
        float bond_r0=nanf("");

        for(const auto &pt : state->polymer_types){
            for(const auto &b : pt.bonds){
                if(isnanf(bond_kappa)){
                    bond_kappa=b.kappa;
                    bond_r0=b.r0;
                }else{
                    if(! (float(b.kappa) == bond_kappa) ){
                        return "All bonds must have same kappa.";
                    }
                    if( !(float(b.r0) == bond_r0) ){
                        return "All bonds must have same kappa.";
                    }
                }
            }

            for(const auto &bp : pt.bond_pairs){
                if(bp.theta0!=0){
                    return "All bond pairs must be straight.";
                }
                if(round(bp.kappa)!=bp.kappa){
                    return "All bond pairs must have integer kappa.";
                }
            }
        }

        return std::string();
    }

    virtual void Attach(WorldState *state)
    {
        m_state=state;
        if(state){
            check_constraints_and_setup();
        }else{
            m_cells.clear();
        }
    }

    virtual void Run(unsigned nSteps) override
    {
        import_beads();

        for(unsigned i=0; i<nSteps; i++){
            step();
            for(const auto &cell : m_cells){
                for(const auto &b : cell.beads){
                    assert( floor(b.x) == cell.location );
                }
            }
        }

        export_beads();
    }

private:
    WorldState *m_state;

    unsigned m_numBeadTypes;
    vec3i_t m_box;

    struct cached_bond_t
    {
        uint32_t bead_hash;
        vec3f_t x;
    };

    struct cell_t
    {
        // Would be encoded in edges. Not true run-time state
        std::vector<unsigned> neighbours;

        vec3i_t location;
        std::vector<bead_resident_t> beads;
        std::vector<cached_bond_t> cached_bonds;

        std::vector<bead_resident_t> outgoing;
    };

    std::vector<cell_t> m_cells;

    float m_dt;
    float m_inv_root_dt;
    uint64_t m_t_hash;

    // We store sqrt(dissipative), to avoid sqrt in inner core maths.
    std::vector<InteractionStrength> m_interactions;

    float m_bond_kappa;
    float m_bond_r0;

    unsigned cell_location_to_cell_index(const vec3i_t &x)
    {
        return x[0] + m_box[0]*x[1] + m_box[0]*m_box[1]*x[2];
    }



    template<class TBeadResident>
    void import_bead(TBeadResident &bb, const Bead &b)
    {
        const PolymerType &pt = m_state->polymer_types.at(b.polymer_type);
        bool is_monomer=pt.bead_types.size()==1;

        bb.id=make_bead_id(is_monomer, b.polymer_id, b.polymer_offset, b.bead_type);
        vec3_copy(bb.x, b.x);
        vec3_copy(bb.v, b.v);
        vec3_copy(bb.f, b.f);
        for(int i=0; i<MAX_BONDS_PER_BEAD; i++){
            bb.bond_partners[i]=-1;
        }
        for(int i=0; i<MAX_ANGLE_BONDS_PER_BEAD; i++){
            bb.angle_bonds[i].partner_head=-1;
            bb.angle_bonds[i].partner_tail=-1;
            bb.angle_bonds[i].kappa=0;
        }

        // Build up the bond info
        if(!is_monomer){
            int nbonds=0;
            // TODO: this is slower than needed. Could be done once.
            for(const auto &bond : pt.bonds){
                if(bond.bead_offset_head==b.polymer_offset){
                    bb.bond_partners[nbonds]=bond.bead_offset_tail;
                    nbonds++;
                }
                if(bond.bead_offset_tail==b.polymer_offset){
                    bb.bond_partners[nbonds]=bond.bead_offset_head;
                    nbonds++;
                }
            }
            require(nbonds <= MAX_BONDS_PER_BEAD, "Too many bonds per bead."); // Should already have been chcked.

            int nangles=0;
            for(const auto &bond_pair : pt.bond_pairs){
                unsigned middle=pt.bonds[bond_pair.bond_offset_head].bead_offset_tail;
                if(b.polymer_offset == middle){
                    bb.angle_bonds[nangles].partner_head=pt.bonds[bond_pair.bond_offset_head].bead_offset_head;
                    bb.angle_bonds[nangles].partner_tail=pt.bonds[bond_pair.bond_offset_tail].bead_offset_tail;
                    bb.angle_bonds[nangles].kappa=(unsigned)bond_pair.kappa;
                    bb.angle_bonds[nangles]._pad_=0;
                    require(bond_pair.theta0==0, "Assume straight bonds."); // Should have been checked earlier
                    nangles++;
                }
            }
            require(nangles <= MAX_ANGLE_BONDS_PER_BEAD, "Too many angles per bead."); // Should already have been chcked.
        }
    }

    void require(bool cond, const char *msg)
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
            m_box[i]=m_state->box[i];
        }

        m_bond_kappa=nanf("");

        // Validate all bonds against assumptions
        for(const auto &pt : m_state->polymer_types){
            std::vector<int> nbonds(pt.bead_types.size(), 0);
            for(const auto &bond : pt.bonds){
                // This could both be changed, with a bit more state and work at run-time.
                if(isnanf(m_bond_kappa)){
                    m_bond_kappa=bond.kappa;
                    m_bond_r0=bond.r0;
                }else{
                    require( float(bond.kappa) == m_bond_kappa, "This method assumes all bonds have the same kappa." );
                    require( float(bond.r0) == m_bond_r0, "This method assumes all bonds have the same r0.");
                }

                nbonds[bond.bead_offset_head]++;
                require(nbonds[bond.bead_offset_head]<=MAX_BONDS_PER_BEAD, "Too many bonds per bead.");
                
                nbonds[bond.bead_offset_tail]++;
                require(nbonds[bond.bead_offset_tail]<=MAX_BONDS_PER_BEAD, "Too many bonds per bead.");
            }

            std::vector<int> nangles(pt.bead_types.size(), 0);
            for(const auto &bond_pair : pt.bond_pairs){
                // This could both be changed, with a bit more state and work at run-time.
                require( bond_pair.theta0==0, "This method assumes all bond pairs are straight.");
                require( bond_pair.kappa<255 && bond_pair.kappa==round(bond_pair.kappa), "This method assumes bond pair kappas are small integers.");
                
                unsigned middle=pt.bonds.at(bond_pair.bond_offset_head).bead_offset_tail;
                nangles[middle]++;
                require(nangles[middle] <= MAX_ANGLE_BONDS_PER_BEAD, "Too many angles per bead.");
            }
        }

        m_numBeadTypes=m_state->bead_types.size();
        m_interactions=m_state->interactions;
        for(auto &ii : m_interactions){
            ii.dissipative=sqrt(ii.dissipative);
        }
        m_dt=m_state->dt;
        m_inv_root_dt=1.0/sqrt(m_state->dt);

        m_cells.resize( m_box[0] * m_box[1] * m_box[2] );
        vec3i_t pos;
        for(pos[0] = 0 ; pos[0] < m_box[0] ; pos[0]++){
            for(pos[1] = 0 ; pos[1] < m_box[1] ; pos[1]++){
                for(pos[2] = 0 ; pos[2] < m_box[2] ; pos[2]++){
                    cell_t &cell = m_cells.at(cell_location_to_cell_index(pos));
                    cell.location=pos;
                    cell.beads.clear();
                    cell.neighbours.clear();
                    cell.outgoing.clear();

                    for_each_point_in_box({-1,-1,-1}, {2,2,2}, [&](const vec3i_t &dir){
                        vec3i_t adj;
                        for(unsigned d=0; d<3; d++){
                            adj[d] = ( (cell.location[d]+dir[d]==m_box[d]) ? -m_box[d] : 0 ) + ( (cell.location[d]+dir[d]==-1) ? +m_box[d] : 0 );
                        }
                        assert(cell.neighbours.size()<27);
                        unsigned neighbour_index=cell_location_to_cell_index(cell.location+dir+adj);
                        cell.neighbours.push_back(neighbour_index);
                    });
                }
            }
        }

        for(const auto &b : m_state->beads){
            bead_resident_t bb;
            import_bead(bb, b);
            vec3i_t cell_location=floor(b.x);
            unsigned cell_index=cell_location_to_cell_index(cell_location);
            m_cells.at(cell_index).beads.push_back(bb);
        }
    }

    void import_beads()
    {
        unsigned done=0;
        for(unsigned cell_index=0; cell_index<m_cells.size(); cell_index++){
            cell_t &c = m_cells.at(cell_index);

            auto it = c.beads.begin();
            while(it!=c.beads.end()){
                auto &bb=*it;
                auto bead_id=m_state->polymers.at(bb.id.get_polymer_id()).bead_ids.at(bb.id.get_polymer_offset());
                const auto &b=m_state->beads.at(bead_id);
                bb.x=vec3f_t(b.x);
                bb.v=vec3f_t(b.v);
                bb.f=vec3f_t(b.f);
                
                unsigned true_index=cell_location_to_cell_index(floor(b.x));
                if(true_index != cell_index){
                    m_cells.at(true_index).beads.push_back(bb);
                    // Count as done if target is an already processed cell, otherwise we'll process again at some point
                    done += (true_index < cell_index); 
                    it=c.beads.erase(it);
                }else{
                    done += 1;
                    ++it;
                }
            }
        }

        if(done!=m_state->beads.size()){
            throw std::runtime_error("World state is corrupted.");
        }
    }

    void export_beads()
    {
        unsigned done=0;
        for(unsigned cell_index=0; cell_index<m_cells.size(); cell_index++){
            const cell_t &c = m_cells.at(cell_index);

            for(const auto &bb : c.beads){
                auto bead_id=m_state->polymers.at(bb.id.get_polymer_id()).bead_ids.at(bb.id.get_polymer_offset());
                auto &b=m_state->beads.at(bead_id);
                b.x=vec3r_t(bb.x);
                b.v=vec3r_t(bb.v);
                b.f=vec3r_t(bb.f);
            }
        }
    }

    vec3i_t world_pos_to_cell_pos(const vec3r_t &pos) const
    { return floor(pos); }

    virtual void step()
    {
        double dt=m_state->dt;

        m_t_hash = next_t_hash(m_state->seed);

        // Move the beads, and then assign to cells based on x(t+dt)
        std::vector<bead_resident_t> moved; // Temporarily hold those that moved
        for(auto &c : m_cells){
            c.cached_bonds.clear();

            unsigned i=0;
            while(i<c.beads.size()){
                auto &b = c.beads[i];
                // Actual position update, clears force
                dpd_maths_core_half_step::update_pos(dt, m_box, b);

                // Check it is still in the right cell
                vec3i_t true_loc=floor(b.x);
                if(true_loc != c.location){
                    // Movement is fairly unlikely
                    moved.push_back(b);
                    if(i+1<c.beads.size()){
                        std::swap(c.beads[i], c.beads.back());
                    }
                    c.beads.resize(c.beads.size()-1);
                }else{
                    ++i;
                }
            }
        }
        // Anything that moved then goes into the right cell
        // This is equivalent to sending messages
        for(const auto &b : moved){
            unsigned true_index=cell_location_to_cell_index(floor(b.x));
            m_cells.at(true_index).beads.push_back(b);
        }

        // Calculate all the non-angle forces
        
        for(cell_t &cell : m_cells){
            assert(cell.neighbours.size()==27);
            unsigned tmp=0;
            for(unsigned neighbour_index : cell.neighbours){
                cell_t &neighbour = m_cells[neighbour_index];
                //std::cerr<<m_state->t<<", i="<<tmp<<", cell="<<cell.location<<", nbor="<<neighbour.location<<"\n";
                tmp++;

                bool self=(&neighbour==&cell);
                // Stay on one neighbour so we can tell if it must be cached
                // We do wrapping on a per bead basis, to simulate receiving each one
                for(unsigned j=0; j<neighbour.beads.size(); j++){
                    bool cache_neighbour=false;
                    const bead_view_t &neighbour_bead=neighbour.beads[j];

                    vec3f_t neighbour_x = neighbour_bead.x;
                    vec3i_t neighbour_cell_pos = floor(neighbour_x);
                    bool wrapped=false;
                    for(int d=0; d<3; d++){
                        if(cell.location[d]==0 && neighbour_cell_pos[d]==m_box[d]-1){
                            neighbour_x[d] -= m_box[d];
                            wrapped=true;
                        }else if(cell.location[d]==m_box[d]-1 && neighbour_cell_pos[d]==0){
                            neighbour_x[d] += m_box[d];
                            wrapped=true;
                        }
                    }
                    if(wrapped){
                        //std::cerr<<"Wrapped: "<<neighbour_bead.x<<" -> "<<neighbour_x<<"\n";
                    }

                    for(unsigned i=0; i<cell.beads.size(); i++){
                        // This implicitly interacts each bead with itself, which is handled with a
                        // distance check in calc_force.
                        cache_neighbour |= calc_force( cell.beads[i], neighbour_bead, neighbour_x);
                    }

                    if(cache_neighbour){
                        uint32_t hash_code=neighbour_bead.id.get_hash_code();
                        cell.cached_bonds.push_back({hash_code, neighbour_x});
                        //std::cerr<<"Caching, t="<<m_state->t<<", cell="<<cell.location<<", "<<neighbour_bead.id.get_polymer_id()<<", poff="<<neighbour_bead.id.get_polymer_offset()<<", hash="<<hash_code<<"\n";
                    }
                }
            }
            for(const auto &b : cell.beads){
                //std::cerr<<"  bpo="<<b.get_hash_code()<<", fdpd="<<b.f<<"\n";
            }
        }

        // At some point we have received and cached all angle bond updates.
        // This must be done pre hardware idle, as any forces need to be send out
        for(cell_t &cell : m_cells){
            std::vector<force_input_t> forces;
            for(auto &b : cell.beads){
                update_bead_angle(cell.cached_bonds, b, forces);
            }
            for(unsigned neighbour_index : cell.neighbours){
                on_recv_forces(m_cells[neighbour_index], forces);
            }
            cell.cached_bonds.clear();
        }

        // Finally do the velocity updates together (hardware idle)
        for(cell_t &cell : m_cells){
            for(auto &b : cell.beads){
                dpd_maths_core_half_step::update_mom(dt, b);
            }
        }

        m_state->t += m_state->dt;
    }

    // Note that b_x has already been adjusted for wrapping
    bool calc_force(bead_resident_t &a, const bead_view_t &b, const vec3f_t &b_x)
    {
        vec3f_t dx=a.x-b_x;
        float dr_sqr=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        if(dr_sqr >=1 || dr_sqr < 1e-5){ // The min threshold avoid large forces, and also skips self-interaction
            return false;
        }
        float dr=sqrt(dr_sqr);

        //std::cerr<<"dut : a="<<a.get_bead_id()<<", b="<<b.get_bead_id()<<", dx="<<dx<<"\n";

        float kappa, r0;
        int angle_bond_count;
        get_bond_info(
            a,
            b.id,
            angle_bond_count,
            kappa,
            r0
        );

        vec3f_t f;
        dpd_maths_core_half_step::calc_force(
            m_inv_root_dt,
            [&](unsigned a, unsigned b){ return m_interactions[m_numBeadTypes*a+b].conservative; },
            // Note that this is teh sqrt of the dissipative
            [&](unsigned a, unsigned b){ return m_interactions[m_numBeadTypes*a+b].dissipative; },
            m_t_hash,
            dx, dr,
            kappa, r0, 
            a, b,
            f
        );
        //std::cerr<<"Dut: t_hash="<<m_t_hash<<", h="<<a.id.get_polymer_id()<<", dx="<<dx<<", dr="<<dr<<", f="<<f<<"\n";

        a.f += f;
        return angle_bond_count > 0;
    }

    // TCache is a vector-like container of cached bonds
    // TOutgoing is a container that supports push_back of a force_input_t
    template<class TCache,class TOutgoing>
    void update_bead_angle(TCache &cache, bead_resident_t &bead, TOutgoing &outgoing)
    {
        auto find_cached_pos=[&](uint32_t target_hash) -> const cached_bond_t *
        {
            for(unsigned i=0; i<cache.size(); i++){
                if(cache[i].bead_hash==target_hash){
                    return &cache[i];
                }
            }
            assert(0); // Should be impossible to get here unless bond has snapped
            return 0;
        };

        for(int i=0; i<MAX_ANGLE_BONDS_PER_BEAD; i++){
            if(bead.angle_bonds[i].partner_head==0xFF){
                break;
            }
            auto head_hash=bead.id.make_hash_from_offset(bead.angle_bonds[i].partner_head);
            auto tail_hash=bead.id.make_hash_from_offset(bead.angle_bonds[i].partner_tail);
            const auto *head=find_cached_pos(head_hash);
            const auto *tail=find_cached_pos(tail_hash);
            
            // The cache copies should already have wrapping applied
            auto first=bead.x - head->x;
            auto second=tail->x - bead.x;

            float FirstLength   = first.l2_norm();
            float SecondLength  = second.l2_norm();
            assert(FirstLength < 1);
            assert(SecondLength < 1);

            vec3f_t headForce, middleForce, tailForce;

            dpd_maths_core_half_step::calc_angle_force(
                (float)bead.angle_bonds[i].kappa, 0.0f, 0.0f,
                first, FirstLength,
                second, SecondLength,
                headForce, middleForce, tailForce
            );

            //std::cerr<<"Dut: fh="<<headForce<<", fm="<<middleForce<<", ft="<<tailForce<<"\n";

            bead.f += middleForce;

            outgoing.push_back( { head_hash, headForce } );
            outgoing.push_back( { tail_hash, tailForce } );
        }
    }

    void on_recv_forces(cell_t &cell, const std::vector<force_input_t> &forces)
    {
        for(auto &b : cell.beads){
            auto bpo = b.id.get_hash_code();
            for(const auto &f : forces){
                if(bpo == f.target_hash){
                    //std::cerr<<"  bpo="<<bpo<<", fdpd="<<b.f<<", fangle="<<f.f<<", f="<<b.f+f.f<<"\n";
                    b.f += f.f;
                }
            }
        }
    }

};

#endif
