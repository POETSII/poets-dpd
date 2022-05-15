#ifndef avx2_dpd_engine_hpp
#define avx2_dpd_engine_hpp

#include <immintrin.h>

#include <atomic>
#include <mutex>
#include <array>
#include <cmath>

#include "dpd/core/dpd_engine.hpp"
#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>

template<unsigned MAX_LOCAL=8>
class AVX2DPDEngine
    : public DPDEngine
{
private:
    static_assert(MAX_LOCAL==8 || MAX_LOCAL==16);

    static const unsigned MAX_BEAD_TYPES = 16;

    static __attribute__((always_inline)) inline void declare(bool cond)
    {
        assert(cond);
        if(!cond){
            __builtin_unreachable();
        }
    }

    // Copies between n and MAX_LOCAL words from src to dst.
    // It is guaranteed that reading MAX_LOCAL values from src is valid,
    // and that writing MAX_LOCAL values is also valid.
    // The reason for all this is:
    // - To provide explicit length hints to the compiler, possibly avoiding normal memcpy.
    // - so that we can explicitly convert into single AVX reads and writes to implement it, if profiling suggests it would be a good idea.
    template<class T>
    static void copy_max_local_words(
        T *dst,
        const T *src,
        unsigned n
    ){
        /*static_assert(sizeof(T)==4);
        declare(n<=MAX_LOCAL);
        std::copy(src, src+n, dst);

        #ifndef NDEBUG
        // In debug we always read and write the extra
        std::copy(src+n, src+MAX_LOCAL, dst+n);
        #endif*/

        static_assert(sizeof(T)*MAX_LOCAL==32 || sizeof(T)*MAX_LOCAL==64);
        __m256i val = _mm256_loadu_si256((const __m256i *)src);
        _mm256_storeu_si256((__m256i*)dst, val);

        static_assert(sizeof(T)==4);
        if(MAX_LOCAL>8 && n>8){
            val = _mm256_loadu_si256((const __m256i *)(src+8));
            _mm256_storeu_si256((__m256i*)(dst+8), val);
        }
    }

    static void copy_and_add_max_local_words(
        float *dst,
        const float *src,
        unsigned n,
        float inc
    ){
        /*static_assert(sizeof(T)==4);
        declare(n<=MAX_LOCAL);
        for(unsigned i=0; i<n; i++){
            dst[i] = src[i] + inc;
        }

        #ifndef NDEBUG
        // In debug we always read and write the extra
        for(unsigned i=n; i<MAX_LOCAL; i++){
            dst[i] = src[i] + inc;
        }
        #endif
        */

        static_assert(sizeof(float)*MAX_LOCAL==32 || sizeof(float)*MAX_LOCAL==64);

        __m256 val = _mm256_loadu_ps(src);
        __m256 xx = _mm256_broadcast_ss(&inc);
        val=_mm256_add_ps(val, xx);
        _mm256_storeu_ps(dst, val);

        if(MAX_LOCAL>8 && n>8){
            val = _mm256_loadu_ps(src+8);
            val=_mm256_add_ps(val, xx);
            _mm256_storeu_ps(dst+8, val);
        }
    }

    struct Cell
    {
        vec3i_t location;
        bool is_on_boundary;
        std::array<uint32_t,27> nhood_info; // Indexes of neighbouring cells. If is_on_boundary is true, then upper 6 bits say what wrapping to apply
        WorldState *world; // Convenience member for debug assertions

        std::atomic<unsigned> local_nhood_count;

        // These are stored in structure of arrays form to make building full nhoods faster
        uint32_t ids[MAX_LOCAL];
        float positions[3][MAX_LOCAL];
        float velocities[3][MAX_LOCAL];        
        float forces[3][MAX_LOCAL];
        uint32_t bead_indices[MAX_LOCAL]; // Store the linear index rather than pointer to bead.

        Cell()
            : local_nhood_count(0)
        {}

        Cell(const Cell &x)
            : location(x.location)
            , is_on_boundary(x.is_on_boundary)
            , nhood_info(x.nhood_info)
            , world(x.world)
            , local_nhood_count(0)
        {
            if(x.local_nhood_count!=0){
                throw std::runtime_error("Only copy empty cells.");
            }
        }

        Cell &operator=(const Cell &o)
        {
            if(o.local_nhood_count!=0){
                throw std::runtime_error("Only copy empty cells.");
            }
            location=o.location;
            is_on_boundary=o.is_on_boundary;
            nhood_info=o.nhood_info;
            local_nhood_count=0;
            world=o.world;
            return *this;
        }

        void nhood_clear()
        {
            local_nhood_count=0;
        }

        template<class T>
        unsigned add_to_local(uint32_t bead_index, uint32_t id, T position, T velocity, T force)
        {
            assert(world->beads[bead_index].get_hash_code().hash==id);

            unsigned index=local_nhood_count.fetch_add(1, std::memory_order::memory_order_relaxed);
            if(index>=MAX_LOCAL){
                throw std::runtime_error("MAX_LOCAL too small.");
            }

            bead_indices[index]=bead_index;
            ids[index]=id;
            for(int i=0; i<3; i++){
                positions[i][index] = position[i];
                velocities[i][index] = velocity[i];
                forces[i][index] = force[i];
            }

            return index;
        }

        unsigned add_to_local(const Cell &c, unsigned local_index)
        {
            assert(world->beads[c.bead_indices[local_index]].get_hash_code().hash==c.ids[local_index]);

            unsigned index=local_nhood_count.fetch_add(1, std::memory_order::memory_order_relaxed);
            if(index>=MAX_LOCAL){
                throw std::runtime_error("MAX_LOCAL too small.");
            }

            bead_indices[index]=c.bead_indices[local_index];
            ids[index]=c.ids[local_index];
            for(int i=0; i<3; i++){
                positions[i][index] = c.positions[i][local_index];
                velocities[i][index] = c.velocities[i][local_index];
                forces[i][index] = c.forces[i][local_index];
            }

            return index;
        }

        unsigned add_to_local(Bead *b)
        {
            return add_to_local(b->bead_id, b->get_hash_code().hash, b->x, b->v, b->f);
        }

        // Write the current state of local out to the owning bead
        void sync_local_to_bead(WorldState *state, unsigned index)
        {
            assert(index < local_nhood_count);

            Bead &b=state->beads[bead_indices[index]];
            assert(b.bead_id==bead_indices[index]);
            assert(b.get_hash_code().hash==ids[index]);
            for(int i=0; i<3; i++){
                b.x[i] = positions[i][index];
                b.v[i] = velocities[i][index];
                b.f[i] = forces[i][index];
            }
        }
    };

    // Used to collect the full surrounding neighbourhood.
    // In structure of arrays form so that we can use AVX to do 8 interactions at once.
    // There is typically a lot of space between each array, but that is not too
    // much of a problem from a cache point of view.
    // This structure is ~(4+3*4+3*4)*MAX_LOCAL*32 = 28*32*MAX_LOCAL=896*MAX_LOCAL bytes
    // With MAX_LOCAL=8, that's about 7200 bytes.
    // So it covers two or three pages, but the entries should stay hot.
    struct FullNHood
    {
        static const unsigned MAX_COUNT = MAX_LOCAL*32;

        unsigned count;

        // Alignment to 64 bytes works for both AVX2 and AVX512
        uint32_t __attribute__ ((aligned (64))) ids[MAX_COUNT];
        float __attribute__ ((aligned (64))) positions[3][MAX_COUNT];
        float __attribute__ ((aligned (64))) velocities[3][MAX_COUNT];

        // This is mutable working area, used to get around the
        // problem of gcc not inferring shuffles when auto vectorising;
        float __attribute__ ((aligned (64))) conStrength[MAX_COUNT];
        float __attribute__ ((aligned (64))) sqrtDissStrength[MAX_COUNT];


        void add(const Cell &o)
        {
            assert(count+o.local_nhood_count <= MAX_COUNT);
            copy_max_local_words(ids+count, o.ids, o.local_nhood_count);
            for(int i=0; i<3; i++){
                copy_max_local_words(positions[i]+count, o.positions[i], o.local_nhood_count);
                copy_max_local_words(velocities[i]+count, o.velocities[i], o.local_nhood_count);
            }
            count += o.local_nhood_count;
        }

        void add_wrapped(const Cell &o, uint32_t bits, float dims[3])
        {
            assert(count+o.local_nhood_count <= MAX_COUNT);
            
            for(int i=0; i<3; i++){
                uint32_t sel=(bits>>(2*i))&0x3;
                if(sel==0){
                    copy_max_local_words(positions[i]+count, o.positions[i], o.local_nhood_count);
                }else{
                    assert( (sel==2) || (sel==1) ); // We have a boundary on one side or the other
                    float dx= (sel==1) ? -dims[i] : +dims[i];
                    
                    copy_and_add_max_local_words(positions[i]+count, o.positions[i], o.local_nhood_count, dx);
                }
            }

            copy_max_local_words(ids+count, o.ids, o.local_nhood_count);
            for(int i=0; i<3; i++){
                copy_max_local_words(velocities[i]+count, o.velocities[i], o.local_nhood_count);
            }
            count += o.local_nhood_count;
        }
    };

    struct polymer_bead_cache
    {
        // Spatial position
        vec3f_t x;
        // Which cell it is currently in
        uint32_t cell_index : 24;
        uint32_t cell_offset : 8;
    };

    std::vector<Cell> m_cells, m_cells_next;
    // This has one location per bead, but only polymers are actually populated
    std::vector<polymer_bead_cache> m_polymer_positions; 
    int m_dims[3];
    float m_dimsf[3];
    float m_dt;
    uint64_t m_t_hash;
    float m_scale_inv_sqrt_dt;
    float m_conStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES];
    float m_sqrtDissStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES];

    WorldState *m_world;

    template<class T3>
    unsigned get_cell_index(T3 pos) const
    {
        return pos[0] + pos[1]*m_dims[0] + pos[2]*m_dims[0]*m_dims[1];
    }

    vec3i_t get_cell_pos(unsigned index) const
    {
        return { int(index % m_dims[0]), int( (index / m_dims[0] ) % m_dims[1] ) , int( (index / (m_dims[0]*m_dims[1])) ) };
    }

    Cell &get_cell_by_index(vec3i_t index)
    {
        return m_cells[get_cell_index(index)];
    }

    using range_t = tbb::blocked_range<unsigned>;

    void __attribute__((noinline)) import_beads_range(const range_t &r)
    {
        for(unsigned i=r.begin(); i<r.end(); i++){
            const Bead &b=m_world->beads[i];
            vec3i_t cell_location=vec3_floor(b.x);
            unsigned cell_index=get_cell_index(cell_location);
            Cell &cell = m_cells[cell_index];

            // Pre-correct backwards by one step            
            auto velocity_half = b.v - b.f * (m_dt*0.5f);

            unsigned cell_offset = cell.add_to_local( b.bead_id, b.get_hash_code().hash, b.x, velocity_half, b.f);

            if(!b.is_monomer){
                for(int d=0; d<3; d++){
                    m_polymer_positions[b.bead_id].x[d] = b.x[d];
                }
                m_polymer_positions[b.bead_id].cell_index=cell_index;
                m_polymer_positions[b.bead_id].cell_offset=cell_offset;
            }
        }
    }

    void import_beads()
    {
        tbb::parallel_for<range_t>(range_t(0, m_world->beads.size(), 1024), [&](const range_t &r){
            import_beads_range(r);
        });
    }

    void __attribute__((noinline)) export_beads_range(const range_t &r)
    {
        for(unsigned ci=r.begin(); ci<r.end(); ci++){
            Cell &c=m_cells[ci];
            for(unsigned i=0; i<c.local_nhood_count; i++){
                c.sync_local_to_bead(m_world, i);

                // Post-correct velocity forwards, as it hasn't been done
                auto &b=m_world->beads[c.bead_indices[i]];
                assert(!isnanf(b.v[0]));
                b.v = b.v + b.f * (m_dt*0.5f);

                assert(!isnanf(b.v[0]));

            }

            // clear the cell out for the next round
            c.local_nhood_count=0;
        }
    }

    void export_beads()
    {
        tbb::parallel_for(range_t(0, m_cells.size(), 256), [&](const range_t &r){
            export_beads_range(r);
        });
    }

    static void __attribute__((noinline)) calculate_cell_forces(
        float dt,
        uint64_t t_hash,
        float scale_inv_sqrt_dt,
        const float conStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        const float sqrtDissStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        FullNHood &nhood,
        Cell &cell
    ){
        unsigned n=nhood.count;

        for(unsigned local_index=0; local_index<cell.local_nhood_count; local_index++){
            float position[3], velocity[3], force[3];
            for(unsigned d=0; d<3; d++){
                position[d]=cell.positions[d][local_index];
                velocity[d]=cell.velocities[d][local_index];
                force[d]=0.0f;
            }
            uint32_t id=cell.ids[local_index];
            unsigned bead_type=BeadHash{id}.get_bead_type();

            // We have to do this out of the main loop, as it requires things
            // that gcc can't inferr (shuffle)
            const float *conStrengthRow = conStrengthMatrix+bead_type*MAX_BEAD_TYPES;
            const float *sqrtDissStrengthRow = sqrtDissStrengthMatrix+bead_type*MAX_BEAD_TYPES;
            if(MAX_LOCAL==8 && MAX_BEAD_TYPES <= 8){
                __m256 conStrengthVals = _mm256_loadu_ps(conStrengthRow);
                __m256 sqrtDissStrengthVals = _mm256_loadu_ps(sqrtDissStrengthRow);

                unsigned nUp=(n+7)/8;
                for(unsigned i=0; i<nUp; i+=8){
                    static_assert(BeadHash::BEAD_TYPE_BITS==4);
                    __m256i ids=_mm256_loadu_si256((const __m256i*)nhood.ids+i);
                    ids = ids>>28;

                    _mm256_storeu_ps( nhood.conStrength+i, _mm256_permutevar8x32_ps(conStrengthVals, ids));
                    _mm256_storeu_ps( nhood.sqrtDissStrength+i, _mm256_permutevar8x32_ps(sqrtDissStrengthVals, ids));
                }
            }else{
                for(unsigned i=0; i<n; i++){
                    unsigned o_bead_type=BeadHash{nhood.ids[i]}.get_bead_type();
                    nhood.conStrength[i] = conStrengthRow[o_bead_type];
                    nhood.sqrtDissStrength[i] = sqrtDissStrengthRow[o_bead_type];
                }
            }

            for(unsigned i=0; i<n; i++){
                float dx[3];
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    dx[d]=  position[d] - nhood.positions[d][i];
                }
                
                float r2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

                // This is the early stopping point, but it is very likely that
                // at least one lane is active in 8-way SIMD as we have beads
                // from multiple neighbours
                
                float r = sqrtf(r2);
                bool active = (r2<1.0f && r2>0.000001f);
                float inv_r= active ? 1.0f/r : 0;  // Note: this will end up as NaN if r==0, which propagates
                float wr = 1.0f - r;

                float dv[3], ddx[3];
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    ddx[d] = dx[d];
                    dx[d] *= inv_r;
                    dv[d] = velocity[d] - nhood.velocities[d][i];
                }
                
                uint32_t oid=nhood.ids[i];
                unsigned matrix_offset = BeadHash{oid}.get_bead_type();

                float conStrength = nhood.conStrength[i];
                float sqrtDissStrength = nhood.sqrtDissStrength[i];
                
                float conForce = conStrength*wr;
                
                float rdotv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
                float sqrt_gammap = sqrtDissStrength*wr;

                float dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
                float u = dpd_maths_core_half_step::default_hash(t_hash, id, oid);
                float randScale = sqrt_gammap * scale_inv_sqrt_dt;
                float randForce = randScale * u;

                float scaled_force = active ? (conForce + dissForce + randForce) : 0;

                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    force[d] += scaled_force*dx[d];
                }
                assert(!isnanf(force[0]));

                /*if(ForceLogging::logger()){    
                    fprintf(stderr, "Logging\n");        
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dx", 3,ddx);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dr", 1,&r);
                    double tt=scale_inv_sqrt_dt;
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-invrootdt", 1, &tt);
                    double gammap=sqrt_gammap*sqrt_gammap;
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-gammap", 1, &gammap);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-rng", 1, &u);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-con",1, &conForce);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-diss", 1,&dissForce);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-rng-scale",1, &randScale);
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"dpd-rand",1, &randForce);
                    double dpd_force=conForce + dissForce + randForce;
                    double ff[3]={(double)dx[0]*dpd_force,(double)dx[1]*dpd_force,(double)dx[2]*dpd_force};
                    ForceLogging::logger()->LogBeadPairProperty( BeadHash{id},BeadHash{oid},"f_next_dpd", 3,ff);
                }*/
            }

            #pragma GCC unroll 3
            for(int d=0; d<3; d++){
                cell.forces[d][local_index] += force[d];
            }
        }
    }

    #if 0
    static void __attribute__((noinline)) calculate_cell_forces_v2(
        float dt,
        uint64_t t_hash,
        float scale_inv_sqrt_dt,
        const float conStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        const float sqrtDissStrengthMatrix[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
        FullNHood &nhood,
        Cell &cell
    ){
        unsigned n=nhood.count;
        unsigned n_avx2=(n+7)/8;

        for(unsigned local_index=0; local_index<cell.local_nhood_count; local_index++){
            float position[3], velocity[3], force[3];
            for(unsigned d=0; d<3; d++){
                position[d]=cell.positions[d][local_index];
                velocity[d]=cell.velocities[d][local_index];
                force[d]=0.0f;
            }
            uint32_t id=cell.ids[local_index];
            unsigned bead_type=BeadHash{id}.get_bead_type();

            __m256 conStrength = _mm256_loadu_ps(conStrengthMatrix+bead_type*MAX_BEAD_TYPES);
            __m256 sqrtDissStrength = _mm256_loadu_ps(sqrtDissStrengthMatrix+bead_type*MAX_BEAD_TYPES);
            
            for(unsigned i=0; i<n; i+=8){
                static_assert(BeadHash::BEAD_TYPE_BITS==4);
                __m256i ids=_mm256_loadu_si256((const __m256i*)(nhood.ids+i));
                ids = ids>>28;

                __m256 conStrengthSel =  _mm256_permutevar8x32_ps(conStrengthRow, ids);
                __m256 sqrtDissStrengthSel = _mm256_permutevar8x32_ps(sqrtDissStrengthRow, ids);
                
                __m256 dx[3];
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    dx[d] =  position[d] - _mm256_loadu_ps(nhood.positions[d]+i);
                }
                
                __m256 r2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

                // This is the early stopping point, but it is very likely that
                // at least one lane is active in 8-way SIMD as we have beads
                // from multiple neighbours
                
                auto r = _mm256_sqrt_ps(r2);
                auto active_mask = _mm256_and_ps( _mm256_cmp_ps(r2,_mm256_set1_ps(1.0f), _CMP_EQ_OQ) , _mm256_cmp_ps(_mm256_set1_ps(0.000001f), r2, _CMP_EQ_OQ));
                auto inv_r = _mm256_and_ps(active_mask, 1.0f / r);  // Note: this will end up as NaN if r==0, which propagates
                auto wr = 1.0f - r;

                __m256 dv[3];
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    dx[d] *= inv_r;
                    dv[d] = velocity[d] - _mm256_loadu_ps(nhood.velocities[d]+i);
                }
                
                auto oid=_mm256_loadu_si256((__m256i*)(nhood.ids+i));
                
                auto conForce = conStrength*wr;
                
                auto rdotv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
                auto sqrt_gammap = sqrtDissStrength*wr;

                auto dissForce = -sqrt_gammap*sqrt_gammap*rdotv;
                auto u = dpd_maths_core_half_step::default_hash(t_hash, id, oid);
                float randScale = sqrt_gammap * scale_inv_sqrt_dt;
                float randForce = randScale * u;

                float scaled_force = active ? (conForce + dissForce + randForce) : 0;

                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    force[d] += scaled_force*dx[d];
                }
                assert(!isnanf(force[0]));
            }

            #pragma GCC unroll 3
            for(int d=0; d<3; d++){
                cell.forces[d][local_index] += force[d];
            }
        }
    }
    #endif

    void __attribute__((noinline)) move_beads_range(const range_t &r)
    {
        float dt=m_dt;  

        for(unsigned index=r.begin(); index<r.end(); index++){
            auto &cell=m_cells[index];              // Exclusive access
            auto &cell_next=m_cells_next[index];    // Shared access via mutex

            unsigned n=cell.local_nhood_count;
            declare(n <= MAX_LOCAL);
            for(unsigned i=0; i<n; i++){
                float position[3];
                float velocity[3];
                float force[3];
                uint32_t bead_index;
                uint32_t id;

                id=cell.ids[i];
                bead_index=cell.bead_indices[i];
                for(int d=0; d<3; d++){
                    position[d]=cell.positions[d][i];
                    velocity[d]=cell.velocities[d][i];
                    force[d]=cell.forces[d][i];
                }

                // Mom step to get from v(t-dt/2) to v(t) (needs to be pre-corrected before first step)
                for(int d=0; d<3; d++){
                    velocity[d] += force[d] * (dt*0.5f);
                }

                // x(t+dt) = x(t) + v(t)*dt + f(t)*dt*dt/2
                // v(t+dt/2) = v(t) + dt*f(t)/2
                int cell_index[3];
                bool moved=false;
                for(int d=0; d<3; d++){
                    position[d] += velocity[d]*dt + force[d] * dt*(dt*0.5f);
                    velocity[d] += force[d] * (dt*0.5f);
                    force[d] = 0;
                    cell_index[d] = floorf(position[d]);
                    moved |= cell_index[d] != cell.location[d]; 
                }

                Cell *dst=&cell_next;
                uint32_t dst_index=index;
                if(moved){
                    for(int d=0; d<3; d++){
                        if(cell_index[d] < 0){
                            position[d] += m_dimsf[d];
                            cell_index[d] += m_dims[d];
                        }else if(cell_index[d] >= m_dims[d]){
                            position[d] -= m_dimsf[d];
                            cell_index[d] -= m_dims[d];
                        }
                    }
                    dst_index=get_cell_index(cell_index);
                    dst=&m_cells_next[dst_index];
                }

                unsigned cell_offset = dst->add_to_local(bead_index, id, position, velocity, force);

                if(!BeadHash{id}.is_monomer()){
                    for(int d=0; d<3; d++){
                        m_polymer_positions[bead_index].x[d] = position[d];
                    }
                    m_polymer_positions[bead_index].cell_index=dst_index;
                    m_polymer_positions[bead_index].cell_offset=cell_offset;
                }
            }

            cell.local_nhood_count=0;
        }
    }

    vec3f_t calc_distance_from_to(const vec3f_t &base, const vec3f_t &other) const
    {
        vec3f_t res;
        for(unsigned i=0; i<3; i++){
            float len=m_dimsf[i];
            float dx=other[i]-base[i];
            int wrap_neg = dx > len*0.5;
            int wrap_pos = dx < -len*0.5;
            res[i] = dx + (wrap_pos - wrap_neg) * len;
        }
        return res;
    }

    void __attribute__((noinline)) update_polymer_range(const range_t &r)
    {
        for(unsigned i=r.begin(); i<r.end(); i++){
            const auto &polymer=m_world->polymers[i];
            const auto &pt = m_world->polymer_types[polymer.polymer_type];
            const auto &bead_ids=polymer.bead_ids; 

            // Hookean bonds
            for(const auto &b : pt.bonds){
                auto head_index = bead_ids[b.bead_offset_head];
                auto tail_index = bead_ids[b.bead_offset_tail];

                auto head_pos=m_polymer_positions[head_index];
                auto tail_pos=m_polymer_positions[tail_index];

                vec3f_t dx;
                for(int d=0; d<3; d++){
                    dx[d] = tail_pos.x[d] - head_pos.x[d];
                    if(dx[d] > 2.5){
                        dx[d] -= m_dims[d];
                    }else if(dx[d] < -2.5){
                        dx[d] += m_dims[d];
                    }
                }               
                float r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                float r=sqrtf(r2);
                float inv_r=1.0f/r;
                vec3f_t force_head;
                dpd_maths_core::calc_hookean_force(
                    float(b.kappa), float(b.r0), dx, r, inv_r, force_head
                );

                Cell &head_cell=m_cells[head_pos.cell_index];
                Cell &tail_cell=m_cells[tail_pos.cell_index];

                for(int d=0; d<3; d++){
                    head_cell.forces[d][head_pos.cell_offset] += force_head[d];
                    tail_cell.forces[d][tail_pos.cell_offset] -= force_head[d];
                }
            }

            // Angle bonds
            for(const auto &bp : pt.bond_pairs){
                const auto &head_bond=pt.bonds[bp.bond_offset_head];
                const auto &tail_bond=pt.bonds[bp.bond_offset_tail];

                auto head_index = bead_ids[head_bond.bead_offset_head];
                auto centre_index = bead_ids[head_bond.bead_offset_tail];
                auto tail_index = bead_ids[tail_bond.bead_offset_tail];

                auto head_pos=m_polymer_positions[head_index];
                auto centre_pos=m_polymer_positions[centre_index];
                auto tail_pos=m_polymer_positions[tail_index];

                auto first=calc_distance_from_to(head_pos.x, centre_pos.x);
                auto second=calc_distance_from_to(centre_pos.x, tail_pos.x);

                float FirstLength   = first.l2_norm();
                float SecondLength  = second.l2_norm();

                vec3f_t headForce, middleForce, tailForce;

                dpd_maths_core_half_step::calc_angle_force(
                    (float)bp.kappa, cosf(bp.theta0), sinf(bp.theta0), 
                    first, FirstLength,
                    second, SecondLength,
                    headForce, middleForce, tailForce
                );

                Cell &head_cell=m_cells[head_pos.cell_index];
                Cell &centre_cell=m_cells[centre_pos.cell_index];
                Cell &tail_cell=m_cells[tail_pos.cell_index];

                for(int d=0; d<3; d++){
                    head_cell.forces[d][head_pos.cell_offset] += headForce[d];
                    centre_cell.forces[d][centre_pos.cell_offset] += middleForce[d];
                    tail_cell.forces[d][tail_pos.cell_offset] += tailForce[d];
                }
            }
        }
    }

    void __attribute__((noinline)) update_cell_force_range(const range_t &r){
        // Within this task we have the following ownerships:
        // cell : local nhood : read-only, shared
        // cell_next : local nhood : write-only, shared, access managed by atomic counter

        auto pWorking = std::make_shared<FullNHood>();
        auto &working=*pWorking;

        for(unsigned index=r.begin(); index<r.end(); index++){
            auto &cell=m_cells[index];

            working.count=0;
            if(!cell.is_on_boundary){
                for(unsigned i=0; i<27; i++){
                    // Upper bits of nhood_info will be 0
                    working.add( m_cells[cell.nhood_info[i]] );
                }
            }else{
                for(unsigned i=0; i<27; i++){
                    uint32_t info=cell.nhood_info[i];
                    uint32_t bits=info>>(32-6);
                    if( bits == 0 ){
                        working.add( m_cells[info] );
                    }else{
                        unsigned index=(info << 6) >> 6;
                        working.add_wrapped( m_cells[index], bits, m_dimsf );
                    }
                }
            }
            calculate_cell_forces(m_dt, m_t_hash, m_scale_inv_sqrt_dt, m_conStrengthMatrix, m_sqrtDissStrengthMatrix, working, cell);
        }
    }

    void step_cells()
    {

        using range_t = tbb::blocked_range<unsigned>;

        /*if(ForceLogging::logger()){
            ForceLogging::logger()->SetTime(m_world->t);
            ForceLogging::logger()->LogProperty("dt", 1, &m_world->dt);
            double seed_low=m_world->seed &0xFFFFFFFFul;
            double seed_high=m_world->seed>>32;
            ForceLogging::logger()->LogProperty("seed_lo", 1, &seed_low);
            ForceLogging::logger()->LogProperty("seed_high", 1, &seed_high);
            double t_hash_low=m_t_hash&0xFFFFFFFFul;
            double t_hash_high=m_t_hash>>32;
            ForceLogging::logger()->LogProperty("t_hash_lo", 1, &t_hash_low);
            ForceLogging::logger()->LogProperty("t_hash_high", 1, &t_hash_high);
            for(auto &b : m_world->beads){
                double h=b.get_hash_code().hash;
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(), "b_hash", 1, &h);
                double x[3]={b.x[0],b.x[1],b.x[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"x",3,x);
                double v[3]={b.v[0],b.v[1],b.v[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"v",3,v);
                double f[3]={b.f[0],b.f[1],b.f[2]};
                ForceLogging::logger()->LogBeadProperty(b.get_hash_code(),"f",3,f);
            }
        }*/

        tbb::parallel_for<range_t>( range_t(0, m_cells.size(), 1024), [&](const range_t &r){
            move_beads_range(r);
        });

        std::swap(m_cells, m_cells_next);

        tbb::parallel_for<range_t>( range_t(0, m_world->polymers.size(), 256), [&](const range_t &r){
            update_polymer_range(r);
        });

        // We want a reasonable grain size due to the allocation for working set and
        // the caching opportunities between adjacent cells. Maybe this should be a 3d
        // range?
        tbb::parallel_for<range_t>( range_t(0, m_cells.size(), 256), [&](const range_t &r){
            update_cell_force_range(r);
        }); 
    }
public:

    std::string CanSupport(const WorldState *world) const override
    {
        if(world->bead_types.size() > MAX_BEAD_TYPES){
            return "Too many bead types.";
        }

        for(const auto &pt : world->polymer_types){
            if(pt.bond_pairs.size()){
                return "BondPairs are unreliable for this engine. Need debugging.";
            }
        }

        return {};
    }

    void Attach(WorldState *state) override
    {
        if(!state){
            m_world=nullptr;
            return;
        }

        m_world=state;
        
        m_dt=m_world->dt;
        m_scale_inv_sqrt_dt=pow_half( dpd_maths_core::kT * 24 / m_world->dt);

        // One entry per cell, but only polymers are populated
        m_polymer_positions.resize(state->beads.size());

        state->box.extract(m_dims);
        state->box.extract(m_dimsf);
        for(unsigned i=0; i<state->bead_types.size(); i++){
            for(unsigned j=0; j<state->bead_types.size(); j++){
                m_conStrengthMatrix[i*MAX_BEAD_TYPES+j] = state->interactions[i*state->bead_types.size()+j].conservative;
                m_sqrtDissStrengthMatrix[i*MAX_BEAD_TYPES+j] = sqrtf(state->interactions[i*state->bead_types.size()+j].dissipative);
            }
        }
        unsigned ncells=m_dims[0]*m_dims[1]*m_dims[2];
        m_cells.resize(ncells);
        m_cells_next.resize(ncells);

        tbb::parallel_for<unsigned>(0, ncells, [&](unsigned i){
            vec3i_t location=get_cell_pos(i);
            Cell &c=m_cells[i];
            c.location=location;
            c.is_on_boundary=false;
            for(int d=0; d<3; d++){
                c.is_on_boundary |= (location[d]==0) || (location[d]==m_dims[d]-1);
            }
            c.world=m_world;
            c.local_nhood_count=0;
            int dst=0;
            for(int dx=-1; dx<=+1; dx++){
                int ox=(location[0]+m_dims[0]+dx) % m_dims[0];
                unsigned xbits=(location[0]==0 && dx==-1) ? 1 : (location[0]==m_dims[0]-1 && dx==+1) ? 2 : 0; 
                for(int dy=-1; dy<=+1; dy++){
                    int oy=(location[1]+m_dims[1]+dy) % m_dims[1];
                    unsigned ybits=(location[1]==0 && dy==-1) ? 1 : (location[1]==m_dims[1]-1 && dy==+1) ? 2 : 0; 

                    for(int dz=-1; dz<=+1; dz++){
                        int oz=(location[2]+m_dims[2]+dz) % m_dims[2];
                        unsigned zbits=(location[2]==0 && dz==-1) ? 1 : (location[2]==m_dims[2]-1 && dz==+1) ? 2 : 0; 

                        unsigned bits=(zbits<<4) | (ybits<<2) | xbits;
                        assert( (bits!=0) ? c.is_on_boundary : true);

                        unsigned index=get_cell_index(vec3i_t{ox,oy,oz});
                        c.nhood_info[dst++] = (bits<<(32-6)) | index;
                    }
                }
            }

            m_cells_next[i] = c;
        });
    }

    void Run(unsigned nSteps) override
    {
        assert(m_world);

        import_beads();

        for(unsigned i=0; i<nSteps; i++){
            m_t_hash = get_t_hash(m_world->t, m_world->seed);
            step_cells();
            m_world->t++;
        }

        export_beads();
    }

};

#endif
