#ifndef tbb_dpd_engine_hpp
#define tbb_dpd_engine_hpp

#error "Not complete."

#include "dpd/core/dpd_engine.hpp"

#include "dpd_maths_core_half_step.hpp"

#include "make_nhood.hpp"

#include "vec3.hpp"
#include "hash.hpp"

#include <cassert>
#include <cmath>
#include <array>

#include "tbb/parallel_for.h"

/*
    This method must have lambda=0.5

    This is based on naive_dpd_half_step, but adds the following:
    - Forces are only calculated once, then propagated to neighbours
    - Updates are calculated in parallel using parallel_for

    We use a simple method of finding a factor of each dimension, of
    at least 3. This means that there are no conflicts in each dimension.
    If the factors are f0, f1, f2, then the number of groups will be
    f0*f1*f2. Taking the example of 100^3, we would choose factor 4,
    producing 64 groups, each containing 15625 cells.

    If the user chooses prime-length dimensions, this method breaks down
    catastrophically. However, for binary powers or decimal powers it
    works well.

    Angle bonds are calculated in the centre, then propagated in both directions,
    which is why we need the clear separation.

    Beads are represented using two structures:
    - A full double-precision structure at a fixed location
    - A single-precision shadow containing just id, position, and velocity in the cell
    Most access for DPD forces are through the shadows, with good locality.
    Bonds can be done through the full precision structure.

    Bonds are represented as data-structures moved with the bead, though
    the bond information is in a seperate structure. Both hookeand and
    angle bonds are represented using the same structure. Only one bead
    in a hookean bond has the bond, and it will propagate that force
    over to the other bond.
    
    There is a global index which always points to the current (memory) location of
    the bead, and is updated when they migrate.
*/
class NaiveDPDEngineHalfStep
    : public DPDEngine
{
public:
    virtual void Attach(WorldState *state)
    {
        m_state=state;
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

    double m_inv_root_dt;
    uint64_t m_t_hash;

    const int MAX_BONDS_PER_BEAD=4;
    
    struct bead_full
    {
        uint32_t hash_code;
        uint32_t bead_type;
        vec3r_t x, v, f;

        // All bonds are the same length and strength
        struct{
            int32_t partner1;
            int32_t partner2; 
            float kappa;
            float sin_theta0;
            float sin_theta1;
        }bond_partners;

        uint32_t get_hash_code() const
        { return hash_code; }

        unsigned get_bead_type() const
        { return bead_type &0xF; }
    };
    

    struct cell
    {
        vec3i_t location;
        std::vector<bead_full*> beads;
        std::vector<cell*> forward_neighbours;
    };

    std::vector<cell> m_cells;
    std::vector<bead_full> m_beads; 
    std::vector<std::vector<cell*>> m_conflict_groups;

    static std::vector<std::vector<vec3i_t>> create_conflict_groups(vec3i_t dims)
    {
        vec3i_t mods;
        for(int i=0; i<3; i++){
            mods[i]=3;
            while(dims[i] % mods[i]){
                mods[i]++;
            }
        }

        //pos=x+y*dim_x+z*dim_x*dim_y;
        vec3i_t scales{1, mods[0], mods[0]*mods[1]};

        std::vector<std::vector<vec3i_t>> res;
        res.resize(mods[0]*mods[1]*mods[2]);
        for_each_point_in_box({0,0,0}, dims, [&](const vec3i_t &x){
            vec3i_t group_pos=x % mods;
            unsigned group_index=dot(group_pos, scales);
            res.at(group_index).push_back(x);
        });
        return res;
    }

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
        m_inv_root_dt=recip_pow_half(m_state->dt);
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
        m_t_hash=get_t_hash(m_state->t, m_state->seed);

        for(const auto &group : m_conflict_groups)
        {
            tbb::parallel_for( (size_t)0, group.size(), [&](unsigned i) {
                process_cell(group[i]);
            });
        }
    }

    void process_cell(cell *c)
    {
        process_intra(c->beads);
        for(auto nboor : c->forward_neighbours){

        }
    }

    void process_intra(cell *c)
    {
        for(int i=0; i<c->beads.size()-1; i++){
            for(int j=i+1; j<c->beads.size(); j++){
                apply_dpd_and_hookean_forces(c->beads[i], c->beads[j], c->beads[j]->x);
            }
        }
    }

    void process_inter(cell *a, cell *b)
    {
        vec3i_t grid_dx=b->location-a->location;
        vec3r_t grid_adj;
        for(int i=0; i<3; i++){
            // Either grid_dx[i]==-1 or it is large and negative
            grid_adj[i] = (grid_dx[i] < 0) ? m_dims[i] : 0.0;
        }

        for(int j=0; j<b->beads.size(); j++){
            vec3r_t b_x=b->beads[j]->x + grid_adj

            for(int i=0; i<a->beads.size(); i++){
                apply_dpd_and_hookean_forces(a->beads[i], b->beads[j], b_x);
            }
        }
    }

    void apply_dpd_and_hookean_forces(bead_full *a, bead_full *b, const vec3r_t &b_x)
    {
        vec3r_t dx=a->x-b_x;

        double kappa=0;
        double r0=0.5;

        for(int i=0; i< MAX_BONDS_PER_BEAD; i++){

        }
    }
};

#endif
