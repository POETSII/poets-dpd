#ifndef DPD_engine_avx2_half_merge_hpp
#define DPD_engine_avx2_half_merge_hpp

#include <immintrin.h>

#include <random>

#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/edge_wrap.hpp"
#include "dpd/core/make_nhood.hpp"

#include "dpd/maths/dpd_maths_core_simd.hpp"
#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"
#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include <tbb/parallel_for_each.h>

class DPDEngineAVX2HalfMergeTBB
    : public DPDEngine
{
public:
    static const unsigned MAX_BEADS_PER_CELL = 16;
    static const unsigned MAX_BEAD_TYPES = 8;


     std::string CanSupport(const WorldState *state) const override
    {
        for(int d=0;d<3;d++){
            int l=round(state->box[d]);
            if(l!=state->box[d]){
                return "World bounds must be integer.";
            }
            if( l%4 ){
                return "All dimensions must be a multiple of 4.";
            }
        }

        return DPDEngine::CanSupport(state);
    }

    void Attach(WorldState *s)
    {
        if(!s){
            m_state=0;
        }

        m_state=s;

        m_dims=s->box;
        m_dimsf=s->box;
        m_dt=s->dt;
        m_scale_inv_sqrt_dt = pow_half( dpd_maths_core::kT * 24 / m_state->dt);
        for(unsigned i=0; i<s->bead_types.size(); i++){
            for(unsigned j=0; j<s->bead_types.size(); j++){
                const auto &ii = m_state->interactions[i*m_state->bead_types.size()+j];
                m_conservative_matrix[i*MAX_BEAD_TYPES+j] = (double)ii.conservative;
                m_sqrt_dissipative_matrix[i*MAX_BEAD_TYPES+j] = sqrtf((double)ii.dissipative);
            }
        }

        m_urng.seed(s->seed);

        unsigned num_cells=m_state->box.x[0] * m_state->box.x[1] * m_state->box.x[2];

        auto fwd_rel=make_relative_nhood_forwards(true);
        
        m_cells.resize(num_cells);
        for(unsigned i=0; i<m_cells.size(); i++){
            auto &cell=m_cells[i];

            cell.n=0;
            cell.n_move_pending=0;
            cell.pos=get_cell_pos(i);
            cell.rng_state=make_rng_state(s->seed, i);
            cell.wrap_bits=create_wrap_bits(&m_dims.x[0], &cell.pos.x[0]);

            auto fwd_abs=make_absolute_nhood(fwd_rel, m_dims, cell.pos);
            assert(fwd_abs.size()==13);
            for(unsigned j=0; j<13; j++){
                cell.neighbours[j] = get_cell_index(fwd_abs[j]);
            }
            // Fill in any padding
            for(unsigned j=13; j<std::size(cell.neighbours); j++){
                cell.neighbours[j] = cell.neighbours[13];
            }
        }

        make_conflict_groups();

        import_beads();
    }

    void Run(unsigned nSteps)
    {
        for(unsigned i=0; i<nSteps; i++){
            validate();
            step();
            validate();
        }

        synchronise_beads_out();
        fprintf(stderr, "Max occupancy = %u\n", max_occupancy);
    }


private:
    struct Packed
    {
        uint32_t id;
        float x[3];
        float v[3];
        float f[3];
    };

    struct Cell
    {
        uint32_t n;
        uint32_t n_move_pending; // Beads in [0,n_move_pending) have not yet been moved
        vec3i_t pos;

        Packed beads[MAX_BEADS_PER_CELL];
        __m256i rng_state;
        uint32_t wrap_bits;

        // Only 13 are used. Remainder are padding, but may get fetched by prefetcher
        uint32_t neighbours[16]; 
    };

    vec3i_t m_dims;
    vec3f_t m_dimsf;
    std::vector<Cell> m_cells;
    float m_dt;
    float m_scale_inv_sqrt_dt;
    float m_conservative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
    float m_sqrt_dissipative_matrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES];

    std::mt19937_64 m_urng;

    std::vector<std::vector<std::vector<unsigned>>> m_cell_waves;

    WorldState *m_state=0;

    unsigned max_occupancy=0;

    __m256i make_rng_state(uint64_t seed, unsigned index)
    {
        // HACK! FIX ME
        return _mm256_set_epi32(
            m_urng(), m_urng(), m_urng(), m_urng(),
            m_urng(), m_urng(), m_urng(), m_urng()
        );
    }

    template<class T,class F>
    void parallel_for_each(std::vector<T> &x, unsigned grain, F &&f)
    {
        using range_t=tbb::blocked_range<size_t>;
        tbb::parallel_for(range_t(0,x.size(),grain), [&](const range_t &r){
            for(size_t i=r.begin(); i<r.end(); i++){
                f(x[i]);
            }
        }, tbb::simple_partitioner{});
    }

    template<class F>
    void parallel_for_each_cell_blocked(F &&f)
    {
        for(unsigned i=0; i<m_cell_waves.size(); i++){
            parallel_for_each(m_cell_waves[i], 1, [&](const std::vector<unsigned> &c) {
                for(auto index : c){
                    f(m_cells[index]);
                }
            });
        }
    }

    unsigned get_cell_index(const vec3i_t &pos)
    {
        for(int d=0; d<3; d++){
            assert(0<=pos[d] && pos[d]<m_dims[d] );
        }
        return m_dims[0]*m_dims[1]*pos[2] + m_dims[0] * pos[1] + pos[0];
    }

    vec3i_t get_cell_pos(unsigned index)
    {
        assert(index<m_cells.size());
        vec3i_t res;
        unsigned working=index;
        res[0] = working % m_dims[0];
        working /= m_dims[0];
        res[1] = working % m_dims[1];
        working /= m_dims[1];
        res[2] = working;
        assert(get_cell_index(res)==index);
        return res;
    }

    void validate()
    {
        for(auto &cell : m_cells){
            for(unsigned i=0; i<cell.n; i++){
                auto &bead = cell.beads[i];
                for(int d=0; d<3; d++){
                    assert(!isnanf(bead.f[d]));
                }
            }
        }
    }

    void make_conflict_groups()
    {
        std::unordered_set<unsigned> seen;

        m_cell_waves.clear();
        for(unsigned gx=0; gx<4; gx+=2){
            for(unsigned gy=0; gy<4; gy+=2){
                for(unsigned gz=0; gz<4; gz+=2){
                    // This gives the origin of a 2x2x2 cube within a 4x4x4 block
                    std::vector<std::vector<unsigned>> group;
                    // Loop over all super cells within group
                    for(unsigned ix=gx; ix<(unsigned)m_dims[0]; ix+=4){
                        for(unsigned iy=gy; iy<(unsigned)m_dims[1]; iy+=4){
                            for(unsigned iz=gz; iz<(unsigned)m_dims[2]; iz+=4){
                                std::vector<unsigned> sc;
                                unsigned off=0;
                                for(int lx=0; lx<2; lx++){
                                    for(int ly=0; ly<2; ly++){
                                        for(int lz=0; lz<2; lz++){
                                            unsigned index=get_cell_index({int(ix+lx),int(iy+ly),int(iz+lz)});
                                            if(!seen.insert(index).second){
                                                throw std::runtime_error("Duplicate");
                                            }
                                            sc.push_back(index);
                                        }
                                    }
                                }
                                group.push_back(std::move(sc));
                            }
                        }
                    }
                    m_cell_waves.push_back(std::move(group));
                }
            }
        }
    }
    void import_beads()
    {
        // This one is a bit difficult to parallelise, but it
        // should not be on the critical path for most scenarios.
        for(const auto &b : m_state->beads){
            auto pos=vec3_floor(b.x);
            unsigned dst_index=get_cell_index(pos);
            auto &dst_cell=m_cells[dst_index];

            if(dst_cell.n==MAX_BEADS_PER_CELL){
                throw std::runtime_error("Too many beads in one cell.");
            }

            Packed p;
            p.id=b.get_hash_code().hash;
            b.x.extract(p.x);
            b.v.extract(p.v);
            b.f.extract(p.f);
            for(int d=0; d<3; d++){
                assert(!isnanf(p.f[d]));
            }

            // Pre-correct mom backwards
            dpd_maths_core_half_step_raw::update_mom(-m_dt, p);

            dst_cell.beads[dst_cell.n]=p;
            dst_cell.n_move_pending++;
            dst_cell.n++;
        }
    }

    void synchronise_beads_out()
    {
        parallel_for_each(m_cells, 32, [&](Cell &cell){
            assert(cell.n==cell.n_move_pending);
            for(unsigned i=0; i<cell.n; i++){
                auto &p=cell.beads[i];
                BeadHash bh{p.id};
                const auto &polymer = m_state->polymers.at(bh.get_polymer_id());
                auto bead_index = polymer.bead_ids.at(bh.get_polymer_offset());

                auto &b=m_state->beads[bead_index];
                assert(b.get_hash_code() == bh);

                b.x.assign(p.x);
                b.v.assign(p.v);
                b.f.assign(p.f);

                // Post-correct for missing step
                dpd_maths_core_half_step::update_mom(m_dt, b);
            }
        });
    }

    void convert_AoS_to_SoA_vec8(
        const Cell &in,
        __m256i &bead_type,
        __m256 x[3],
        __m256 v[3]
    )
    {
        assert(in.n <= 8);
        if(in.n > 8){ __builtin_unreachable(); }

        bead_type = _mm256_setzero_si256();
        // Initialise positions somewhere that is guaranteed not to interact.
        // This means force calculations will automatically return 0 for inactive lanes.
        for(int d=0; d<3; d++){
            x[d] = _mm256_set1_ps(-2.0f);
        }

        // Indexing into a _mm256i treats it as 4x64.
        // It might be more efficient to do the floats this way as well.
        uint32_t bead_types_tmp[8]={0,0,0,0, 0,0,0,0 };

        for(unsigned i=0; i<in.n; i++){
            bead_types_tmp[i] = BeadHash{in.beads[i].id}.get_bead_type();
            for(int d=0; d<3; d++){
                assert(!isnanf(in.beads[i].f[d]));
                x[d][i] = in.beads[i].x[d];
                v[d][i] = in.beads[i].v[d];
            }
        }
        bead_type=_mm256_loadu_si256( (const __m256i *)bead_types_tmp );

    }

    virtual void step()
    {
        double dt=m_state->dt;

        for(const auto &c : m_cells){
            assert(c.n==c.n_move_pending);
        }
        validate();

        parallel_for_each_cell_blocked([&](Cell &cell){
            move_cell(cell);
        });

        validate();

        // At this point cell.n==0
        parallel_for_each_cell_blocked([&](Cell &cell){
            step_cell(cell);
        });

        for(const auto &c : m_cells){
            assert(c.n==c.n_move_pending);
        }
        validate();

        m_state->t++;
    }

    void move_cell(Cell &home)
    {
        if(home.wrap_bits){
            move_cell_spec<true>(home);
        }else{
            move_cell_spec<false>(home);
        }
    }

    template<bool IsEdge>
    void move_cell_spec(Cell &home)
    {
        for(int bi=home.n_move_pending-1; bi>=0; bi--){
            auto &b=home.beads[bi];
            //fprintf(stderr, "Bead %u\n", b.id);
            
            dpd_maths_core_half_step_raw::update_mom(m_dt, b);

            // TODO Not sure if specialising here is a good tradeoff. The boundary
            // checks create dependencies, but there are three of them in parallel
            if(IsEdge){
                dpd_maths_core_half_step_raw::update_pos(m_dt, m_dimsf, b);
            }else{
                dpd_maths_core_half_step_raw::update_pos_no_wrap(m_dt, b);
            }

            vec3i_t new_pos;
            vec3_floor(&new_pos[0], &b.x[0]);
            if( new_pos != home.pos ){
                unsigned new_index = get_cell_index(new_pos);

                // This is safe due to waves/groups
                auto &dst=m_cells[new_index];
                if(dst.n == MAX_BEADS_PER_CELL){
                    throw std::runtime_error("Too many beads in one cell.");
                }
                dst.beads[dst.n] = b;
                dst.n += 1;

                home.n -= 1;
                b = home.beads[home.n];             
            }
        }

        home.n_move_pending=0;
    }

    void step_cell(Cell &home)
    {
        if(home.n > 8){
            throw std::runtime_error("Not implemented yet.");
        }

        if(home.n==0){
            return;
        }

        // TODO : remove
        max_occupancy=std::max(max_occupancy, home.n);

        __m256i home_bead_type;
        __m256 home_x[3], home_v[3], home_f[3];
        __m256i rng_state;

        convert_AoS_to_SoA_vec8( home, home_bead_type, home_x, home_v );
        for(int d=0; d<3; d++){
            home_f[d] = _mm256_setzero_ps();
        }
        rng_state=home.rng_state;

        for(int i=0; i<home.n; i++){
            for(int d=0; d<3; d++){
                assert(!isnanf(home_f[d][i]));
                assert(!isnanf(home.beads[i].f[d]));
            }
        }

        step_cell_intra_vec8(home, rng_state, home_bead_type, home_x, home_v, home_f);

        for(int i=0; i<home.n; i++){
            for(int d=0; d<3; d++){
                assert(!isnanf(home_f[d][i]));
                assert(!isnanf(home.beads[i].f[d]));
            }
        }
        
        // flag is mostly true, and there is loads of iteration inside
        // so icache shouldn't be a problem
        if(home.wrap_bits){
            step_cell_inter_vec8<true>(home, rng_state, home_bead_type, home_x, home_v, home_f);
        }else{
            step_cell_inter_vec8<false>(home, rng_state, home_bead_type, home_x, home_v, home_f);
        }

        for(int i=0; i<home.n; i++){
            for(int d=0; d<3; d++){
                assert(!isnanf(home_f[d][i]));
            }
        }
        
        for(unsigned i=0; i<home.n; i++){
            for(int d=0; d<3; d++){
                assert(!isnanf(home.beads[i].f[d]));
                home.beads[i].f[d] += home_f[d][i];
                assert(!isnanf(home.beads[i].f[d]));
            }
        }
        home.rng_state=rng_state;

        // Set it here so that it is ready for the next move round.
        home.n_move_pending = home.n;
    }

    void step_cell_intra_vec8(
        Cell &home,
        __m256i &rng_state,
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],
        __m256 home_f[3]
    ){
        // All home cells are very likely to interact, assuming r=0.5
        // Only if they are in opposite corners to they miss each other.

        assert(home.n<=8);

        __m256i iota=_mm256_set_epi32(7,6,5,4,3,2,1,0);

        for(unsigned j=1; j<home.n; j++){
            auto &other_bead=home.beads[j];

            __m256 f[3];
            if(dpd_maths_core_simd::interact_vec8_to_scalar<MAX_BEAD_TYPES>(
                m_scale_inv_sqrt_dt,
                m_conservative_matrix,
                m_sqrt_dissipative_matrix,
                rng_state,
                
                home_bead_types,
                home_x,
                home_v,

                BeadHash{other_bead.id}.get_bead_type(),
                other_bead.x,
                other_bead.v,

                f
            )){
                // There is no less than in AVX2, so rely on greater than.
                // Create a mask that is true for [0,j) and false for [j,8)
                __m256i active = _mm256_cmpgt_epi32(_mm256_set1_epi32(j),iota);
                for(int d=0; d<3; d++){
                    __m256 forces_active = _mm256_and_ps(f[d], _mm256_castsi256_ps(active));
                    home_f[d] = home_f[d] + forces_active;
                    other_bead.f[d] -= dpd_maths_core_simd::_mm256_reduce_add_ps(forces_active);
                }
            }
        }
    }

    template<bool IsEdgeCell>
    void step_cell_inter_vec8(
        Cell &home,
        __m256i &rng_state,
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],
        __m256 home_f[3]
    ){
        for(unsigned i=0; i<13; i++){
            auto &other_cell=m_cells[home.neighbours[i]];

            for(unsigned j=0; j<other_cell.n; j++){
                auto &other_bead=other_cell.beads[j];

                __m256 f[3]; // Only set if there is an interaction

                float other_bead_x_local[3];
                if(IsEdgeCell){
                    vec3_copy(other_bead_x_local, other_bead.x);
                    do_neighbour_wrap(other_bead_x_local, home.wrap_bits, &m_dimsf.x[0]);
                }
                
                if(dpd_maths_core_simd::interact_vec8_to_scalar<MAX_BEAD_TYPES>(
                    m_scale_inv_sqrt_dt,
                    m_conservative_matrix,
                    m_sqrt_dissipative_matrix,
                    rng_state,
                    
                    home_bead_types,
                    home_x,
                    home_v,

                    BeadHash{other_bead.id}.get_bead_type(),
                    IsEdgeCell ? other_bead_x_local : other_bead.x,
                    other_bead.v,

                    f
                )){
                    for(int d=0; d<3; d++){
                        for(int i=0; i<home.n; i++){
                            assert(!isnanf(f[d][i]));
                        }
                        home_f[d] = home_f[d] + f[d];
                        assert(!isnanf(other_bead.f[d]));
                        other_bead.f[d] -= dpd_maths_core_simd::_mm256_reduce_add_ps(f[d]);
                        assert(!isnanf(other_bead.f[d]));
                    }
                }
            }
        }
    }



};

#endif
