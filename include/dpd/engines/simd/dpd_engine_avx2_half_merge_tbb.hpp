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

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>

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
                // This could be fixed/reduced
                return "All dimensions must be a multiple of 4.";
            }
        }

        for(const auto &pt : state->polymer_types){
            for(const auto &bt : pt.bond_pairs){
                if(bt.theta0!=0.0){
                    return "Only straight bond angles are supported.";
                }
            }
        }

        return DPDEngine::CanSupport(state);
    }

    double GetMaxBondLength() const override
    { return 2.0; }

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

        m_max_polymer_bonds=0;
        m_max_polymer_length=1;
        for(const auto &pt : s->polymer_types){
            m_max_polymer_bonds=std::max<unsigned>(m_max_polymer_bonds, pt.bonds.size());
            m_max_polymer_length=std::max<unsigned>(m_max_polymer_length, pt.bead_types.size());
        }

        m_non_monomer_bead_locations.clear();
        for(int polymer_id=s->polymers.size()-1; polymer_id>=0; polymer_id--){
            const auto &poly=s->polymers[polymer_id];
            if( poly.bead_ids.size() <= 1){
                continue; // It's a monomer
            }
            if(polymer_id >= m_non_monomer_bead_locations.size()){
                // This should happen just once
                m_non_monomer_bead_locations.resize(polymer_id+1);
            }
            m_non_monomer_polymer_ids.push_back(polymer_id);
            m_non_monomer_bead_locations[polymer_id].reset(new Packed*[poly.bead_ids.size()]);
        }

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

    std::vector<std::unique_ptr<Packed*[]>> m_non_monomer_bead_locations;
    std::vector<uint32_t> m_non_monomer_polymer_ids;
    unsigned m_max_polymer_length;
    unsigned m_max_polymer_bonds;

    WorldState *m_state=0;

    __m256i make_rng_state(uint64_t seed, unsigned index)
    {
        // HACK! FIX ME
        return _mm256_set_epi32(
            m_urng(), m_urng(), m_urng(), m_urng(),
            m_urng(), m_urng(), m_urng(), m_urng()
        );
    }

    static const bool no_par=false;

    template<class T,class F>
    void parallel_for_each(std::vector<T> &x, unsigned grain, F &&f)
    {
        if(no_par){
            for(unsigned i=0; i<x.size(); i++){
                f(x[i]);
            }
        }else{
            using range_t=tbb::blocked_range<size_t>;
            tbb::parallel_for(range_t(0,x.size(),grain), [&](const range_t &r){
                for(size_t i=r.begin(); i<r.end(); i++){
                    f(x[i]);
                }
            }, tbb::simple_partitioner{});
        }
    }

    template<class F>
    void parallel_for_each_cell_blocked(F &&f)
    {
        if(no_par){
            for(unsigned i=0; i<m_cell_waves.size(); i++){
                for(const auto & c : m_cell_waves[i]){
                    for(auto index : c){
                        f(m_cells[index]);
                    }
                };
            }
        }else{
            for(unsigned i=0; i<m_cell_waves.size(); i++){
                parallel_for_each(m_cell_waves[i], 1, [&](const std::vector<unsigned> &c) {
                    for(auto index : c){
                        f(m_cells[index]);
                    }
                });
            }
        }
    }

    template<class F>
    void parallel_for_each_cell_group(F &&f)
    {
        if(no_par){
            for(unsigned i=0; i<m_cell_waves.size(); i++){
                for(const auto &c : m_cell_waves[i]) {
                    f(c);
                };
            }
        }else{
            for(unsigned i=0; i<m_cell_waves.size(); i++){
                parallel_for_each(m_cell_waves[i], 1, [&](const std::vector<unsigned> &c) {
                    f(c);
                });
            }
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

        unsigned cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        // TODO: Do something cleverer here.

        unsigned wx=8, wy=8, wz=8;
        if(wx%8){
            wx=4;
        }
        if(wy%8){
            wy=4;
        }
        if(wz%8){
            wz=4;
        }


        // Number of tasks per wave
        unsigned tasks=( (m_dims[0]/wx) * (m_dims[1]/wy) * (m_dims[2]/wz) ) / 8;
        if(tasks>=4*cpus && (m_dims[1]%16)==0){
            tasks /= 2;
            wy *= 2;
        }
        if(tasks>=4*cpus && (m_dims[2]%16)==0){
            tasks /= 2;
            wz *= 2;
        }
        if(tasks>=4*cpus && (m_dims[0]%16)==0){
            tasks /= 2;
            wx *= 2;
        }
        std::cerr<<"wx="<<wx<<", wy="<<wy<<", wz="<<wz<<", tasks_per_wave="<<tasks<<"\n";

        m_cell_waves.clear();
        for(unsigned gx=0; gx<wx; gx+=(wx/2)){
            for(unsigned gy=0; gy<wy; gy+=wy/2){
                for(unsigned gz=0; gz<wz; gz+=wz/2){
                    // This gives the origin of a 2x2x2 cube within a 4x4x4 block
                    std::vector<std::vector<unsigned>> group;
                    // Loop over all super cells within group
                    for(unsigned ix=gx; ix<(unsigned)m_dims[0]; ix+=wx){
                        for(unsigned iy=gy; iy<(unsigned)m_dims[1]; iy+=wy){
                            for(unsigned iz=gz; iz<(unsigned)m_dims[2]; iz+=wz){
                                std::vector<unsigned> sc;
                                unsigned off=0;
                                for(int lx=0; lx<wx/2; lx++){
                                    for(int ly=0; ly<wy/2; ly++){
                                        for(int lz=0; lz<wz/2; lz++){
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
            if(!b.is_monomer){
                std::unique_ptr<Packed*[]> &loc= m_non_monomer_bead_locations.at(b.polymer_id);
                loc[(unsigned)b.polymer_offset] = dst_cell.beads + dst_cell.n;
            }

            dst_cell.n_move_pending++;
            dst_cell.n++;
        }
    }

    void synchronise_beads_out()
    {
        parallel_for_each(m_cells, 256, [&](Cell &cell){
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
        unsigned base,
        __m256i &bead_type,
        __m256 x[3],
        __m256 v[3]
    )
    {
        bead_type = _mm256_setzero_si256();
        // Initialise positions somewhere that is guaranteed not to interact.
        // This means force calculations will automatically return 0 for inactive lanes.
        for(int d=0; d<3; d++){
            x[d] = _mm256_set1_ps(-2.0f);
        }

        // Indexing into a _mm256i treats it as 4x64.
        // It might be more efficient to do the floats this way as well.
        uint32_t bead_types_tmp[8]={0,0,0,0, 0,0,0,0 };

        unsigned todo=std::min(8u, in.n-base);
        for(unsigned i=0; i<todo; i++){
            bead_types_tmp[i] = BeadHash{in.beads[i+base].id}.get_bead_type();
            for(int d=0; d<3; d++){
                assert(!isnanf(in.beads[i+base].f[d]));
                x[d][i] = in.beads[i+base].x[d];
                v[d][i] = in.beads[i+base].v[d];
            }
        }
        bead_type=_mm256_loadu_si256( (const __m256i *)bead_types_tmp );
    }

    void convert_AoS_to_SoA_dual_vec4(
        const Cell &in,
        __m256i &bead_type,
        __m256 x[3],
        __m256 v[3]
    )
    {
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
                x[d][i] = in.beads[i].x[d];
                v[d][i] = in.beads[i].v[d];
            }
        }
        bead_type=_mm256_loadu_si256( (const __m256i *)bead_types_tmp );

        bead_type=_mm256_insertf128_si256(bead_type, _mm256_castsi256_si128(bead_type), 1);
        for(int d=0; d<3; d++){
            x[d]=_mm256_insertf128_ps(x[d], _mm256_castps256_ps128(x[d]), 1);
            v[d]=_mm256_insertf128_ps(v[d], _mm256_castps256_ps128(v[d]), 1);
        }
    }

    virtual void step()
    {
        double dt=m_state->dt;

        for(const auto &c : m_cells){
            assert(c.n==c.n_move_pending);
        }
        validate();

        parallel_for_each_cell_group( [&](const std::vector<unsigned> &c){
            move_cell_group(c);
        });

        validate();

        // At this point cell.n==0
        parallel_for_each_cell_group( [&](const std::vector<unsigned> &c){
            calc_forces_for_group(c);
        });

        for(const auto &c : m_cells){
            assert(c.n==c.n_move_pending);
        }
        validate();

        m_state->t++;
    }

    __attribute__((noinline)) void move_cell_group(const std::vector<unsigned> &indices)
    {
        for(auto i : indices){
            move_cell(m_cells[i]);
        }
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
                if(!BeadHash{b.id}.is_monomer()){
                    auto &loc= m_non_monomer_bead_locations.at(BeadHash{b.id}.get_polymer_id())[BeadHash{b.id}.get_polymer_offset()];
                    assert(loc == &b);
                    loc = dst.beads+dst.n;
                }
                dst.n += 1;

                home.n -= 1;
                b = home.beads[home.n];             
            }
        }

        home.n_move_pending=0;
    }

    // Making this noinline is good for profiling, but drops performance about 5%
    __attribute__((noinline))
    void calc_forces_for_group(const std::vector<unsigned> &indices)
    {
        for(auto i : indices){
            calc_forces_for_cell(m_cells[i]);
        }
    }


    void calc_forces_for_cell(Cell &home)
    {
        if(home.n==0){
            return;
        }

        __m256i home_bead_type;
        __m256 home_x[3], home_v[3], home_f[3];
        __m256i rng_state;
        
        rng_state=home.rng_state;

        if(home.n <= 4){
            convert_AoS_to_SoA_dual_vec4( home, home_bead_type, home_x, home_v );
            for(int d=0; d<3; d++){
                home_f[d] = _mm256_setzero_ps();
            }

            // This still works for dual vec4
            step_cell_intra_vec8(home, 0, rng_state, home_bead_type, home_x, home_v, home_f);

            // flag is mostly true, and there is loads of iteration inside
            // so icache shouldn't be a problem
            if(home.wrap_bits){
                step_cell_inter_dual_vec4<true>(home, rng_state, home_bead_type, home_x, home_v, home_f);
            }else{
                step_cell_inter_dual_vec4<false>(home, rng_state, home_bead_type, home_x, home_v, home_f);
            }
            
            for(unsigned i=0; i<home.n; i++){
                assert(0<=i && i<home.n);
                for(int d=0; d<3; d++){
                    home.beads[i].f[d] += home_f[d][i] + home_f[d][i+4];
                }
            }
        }else{
            for(unsigned base=0; base<home.n; base+=8){
                unsigned upper=std::min(base+8, home.n);
                assert(upper < home.n );
                assert(base - upper <= 8);

                convert_AoS_to_SoA_vec8( home, base, home_bead_type, home_x, home_v );
                for(int d=0; d<3; d++){
                    home_f[d] = _mm256_setzero_ps();
                }

                step_cell_intra_vec8(home, base, rng_state, home_bead_type, home_x, home_v, home_f);

                // flag is mostly true, and there is loads of iteration inside
                // so icache shouldn't be a problem
                if(home.wrap_bits){
                    step_cell_inter_vec8<true>(home, rng_state, home_bead_type, home_x, home_v, home_f);
                }else{
                    step_cell_inter_vec8<false>(home, rng_state, home_bead_type, home_x, home_v, home_f);
                }

                
                for(unsigned i=base; i<upper; i++){
                    assert(0<=i && i<home.n);
                    for(int d=0; d<3; d++){
                        home.beads[i].f[d] += home_f[d][i-base];
                    }
                }
            }
        }

        home.rng_state=rng_state;

        // Set it here so that it is ready for the next move round.
        home.n_move_pending = home.n;
    }

    void step_cell_intra_vec8(
        Cell &home,
        unsigned base, // Offset of home vector within cell
        __m256i &rng_state,
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],
        __m256 home_f[3]
    ){
        // All home cells are very likely to interact, assuming r=0.5
        // Only if they are in opposite corners to they miss each other.

        assert(home.n<=8);

        __m256i iota=_mm256_set_epi32(7,6,5,4,3,2,1,0) + _mm256_set1_epi32(base);

        for(unsigned j=base+1; j<home.n; j++){
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

            // Prefetch ahead by two neighbours. About 5% perf increase on byron
            _mm_prefetch(&m_cells[home.neighbours[i+2]],_MM_HINT_T1); // This will prefetch past the end

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

    
    template<bool IsEdgeCell>
    void step_cell_inter_dual_vec4(
        Cell &home,
        __m256i &rng_state,
        const __m256i home_bead_types,
        const __m256 home_x[3],
        const __m256 home_v[3],
        __m256 home_f[3]
    ){
        assert(home.n <= 4);

        // Declare before "others" array. Hopefully ends up at bottom of stack
        float other_bead_A_x_local[3];
        float other_bead_A_v_local[3];
        float other_bead_B_x_local[3];
        float other_bead_B_v_local[3];
        __m256 f[3]; // Only set if there is an interaction
        std::pair<float,float> fAB; // Needed to deal with over-run

        // We leave extra padding to make prefetches valid, and
        // also to allow us to always work with pairs
        Packed *others[MAX_BEADS_PER_CELL*13+4]; // This is about 1.6KB

        unsigned num_others=0;
        if(0){
            for(unsigned i=0; i<13; i++){
                auto &other_cell=m_cells[home.neighbours[i]];
                for(unsigned j=0; j<other_cell.n; j++){
                    others[num_others+j] = other_cell.beads+j;
                }
                num_others += other_cell.n;
            }
        }else{
            // Break the dependency between loading the counts and immediately branching on them.
            // Has a decent effect, about 5%
            unsigned other_counts[13];
            for(unsigned i=0; i<13; i++){
                other_counts[i] = m_cells[home.neighbours[i]].n;
            }
            for(unsigned i=0; i<13; i++){
                auto &other_cell=m_cells[home.neighbours[i]];
                for(unsigned j=0; j<other_counts[i]; j++){
                    others[num_others+j] = other_cell.beads+j;
                    // This provides a noticeable speed-up, another 3-5%
                    _mm_prefetch(other_cell.beads+j, _MM_HINT_T1);
                }
                num_others += other_cell.n;
            }
        }
        // Make sure end of array is valid for pre-fetch
        for(unsigned i=0; i<4; i++){
            others[num_others+i] = others[num_others-1];
        }

        // If num_others is odd then we do the last "other" twice
        for(unsigned i=0; i< num_others ; i+=2){
            _mm_prefetch(others[i+2], _MM_HINT_T1); // This will prefetch past the end
            _mm_prefetch(others[i+3], _MM_HINT_T1); // This will prefetch past the end

            Packed &other_bead_A=*others[i];
            Packed &other_bead_B=*others[i+1];
            
            if(IsEdgeCell){
                vec3_copy(other_bead_A_x_local, other_bead_A.x);
                do_neighbour_wrap(other_bead_A_x_local, home.wrap_bits, &m_dimsf.x[0]);

                vec3_copy(other_bead_B_x_local, other_bead_B.x);
                do_neighbour_wrap(other_bead_B_x_local, home.wrap_bits, &m_dimsf.x[0]);
            }
                
            if(dpd_maths_core_simd::interact_dual_vec4_to_dual_scalar<MAX_BEAD_TYPES>(
                m_scale_inv_sqrt_dt,
                m_conservative_matrix,
                m_sqrt_dissipative_matrix,
                rng_state,
                
                home_bead_types,
                home_x,
                home_v,

                BeadHash{other_bead_A.id}.get_bead_type(),
                IsEdgeCell ? other_bead_A_x_local : other_bead_A.x,
                other_bead_A.v,

                BeadHash{other_bead_B.id}.get_bead_type(),
                IsEdgeCell ? other_bead_B_x_local : other_bead_B.x,
                other_bead_B.v,

                f
            )){
                for(int d=0; d<3; d++){
                    for(int i=0; i<home.n; i++){
                        assert(!isnanf(f[d][i]));
                    }
                    home_f[d] = home_f[d] + f[d];
                    fAB = dpd_maths_core_simd::mm256_reduce_add_ps_dual_vec4(f[d]);
                    other_bead_A.f[d] -= fAB.first;
                    // If num_others is odd, and i==num_others-1, this will over-add to the last bead
                    other_bead_B.f[d] -= fAB.second; 
                }
            }
        }

        // If num_others is odd, then reverse it back out
        if(num_others&1){
            for(int d=0; d<3; d++){
                others[num_others-1]->f[d] += fAB.second;
            }
        }
    }

    struct BondWorking
    {
        vec3f_t dx;
        float dr;
    };

    void evaluate_bonds(uint32_t polymer_id, BondWorking *bond_working)
    {
        const PolymerType &pt=m_state->polymer_types[m_state->polymers[polymer_id].polymer_type];
        auto &locations=m_non_monomer_bead_locations.at(polymer_id);

        for(unsigned i=0; i<pt.bonds.size(); i++){
            const Bond &bond=pt.bonds[i];
            auto head_x=locations[bond.bead_offset_head]->x;
            auto tail_x=locations[bond.bead_offset_tail]->x;

            vec3f_t dx;
            float dr2=0;
            for(int d=0; d<3; d++){
                dx[d]=head_x[d]-tail_x[d];
                if(dx[d] < -2 ){
                    dx[d] += m_dims[d];
                }else if(dx[d] > 2){
                    dx[d] -= m_dims[d];
                }
                dr2 += dx[d];
            }
            float dr=pow_half(dr2);
        
            vec3f_t head_f, tail_f;
            dpd_maths_core_half_step::calc_hookean_force<false, float, vec3f_t, vec3f_t>(
                bond.kappa, bond.r0,
                dx, dr,
                head_f, tail_f
            );

            vec3_add(locations[bond.bead_offset_head]->f, &head_f[0]);
            vec3_add(locations[bond.bead_offset_tail]->f, &tail_f[0]);

            bond_working[i].dx=dx;
            bond_working[i].dr=dr;
        }

        for(unsigned i=0; i<pt.bond_pairs.size(); i++){
            const BondPair &bp=pt.bond_pairs[i];
            const Bond &bhead=pt.bonds[bp.bond_offset_head], &btail=pt.bonds[bp.bond_offset_tail];

            assert(bp.theta0==0.0);

            float cos_theta0=1.0f;
            float sin_theta0=0.0f;

            vec3f_t head_f, mid_f, tail_f;
            dpd_maths_core_half_step::calc_angle_force<float, vec3f_t, vec3f_t>(
                bp.kappa,
                cos_theta0,
                sin_theta0,
                bond_working[bp.bond_offset_head].dx, bond_working[bp.bond_offset_head].dr,
                bond_working[bp.bond_offset_tail].dx, bond_working[bp.bond_offset_tail].dr,
                head_f, mid_f, tail_f
            );

            vec3_add(locations[bhead.bead_offset_head]->f, &head_f[0]);
            vec3_add(locations[bhead.bead_offset_tail]->f, &mid_f[0]);
            vec3_add(locations[btail.bead_offset_tail]->f, &tail_f[0]);
        }

    }
};

#endif

