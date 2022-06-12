#ifndef avx2_dpd_engine_v2_hpp
#define avx2_dpd_engine_v2_hpp

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

template<class R, class F>
void ParallelFor(const R &r, F &&f)
{
    if(0){
        f(r);
    }else if(0){
        R rr=R(r.begin(), r.end(), 1);
        tbb::parallel_for<R>(r, f);
    }else{
        tbb::parallel_for<R>(r, f);
    }
}

template<class R, class F>
void NonParallelFor(const R &r, F &&f)
{
    f(r);
}

/*
Changes versus AVX2DPDEngine:
- Uses seperate local (exclusive) and incoming pools (shared) in each cell
- Allows local beads to spill over into incoming if space runs out
- Overly complex movement logic.

Currently this has some kind of correctness bug with angle bonds that
I can't diagnose. Calculations are slightly off, but it is non-determinstic
and trying to debug it makes it go away...

AFAIK it works fine for non-bonded and hookean.
*/
class AVX2DPDEngineV2
    : public DPDEngine
{
private:
    static const unsigned MAX_LOCAL = 8;
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
        if(1){
            static_assert(sizeof(T)==4);
            declare(n<=MAX_LOCAL);
            std::copy(src, src+n, dst);

            #ifndef NDEBUG
            // In debug we always read and write the extra
            std::copy(src+n, src+MAX_LOCAL, dst+n);
            #endif
        }else{
            static_assert(sizeof(T)*MAX_LOCAL==32 || sizeof(T)*MAX_LOCAL==64);
            __m256i val = _mm256_loadu_si256((const __m256i *)src);
            _mm256_storeu_si256((__m256i*)dst, val);
            if(MAX_LOCAL > 8 && n>8){
                val = _mm256_loadu_si256(1+(const __m256i *)src);
                _mm256_storeu_si256(1+(__m256i*)dst, val);        
            }
        }
    }

    static void copy_and_add_max_local_words(
        float *dst,
        const float *src,
        unsigned n,
        float inc
    ){
        if(1){
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
        }else{
            static_assert(sizeof(float)*MAX_LOCAL==32 || sizeof(float)*MAX_LOCAL==64);

            __m256 val = _mm256_loadu_ps(src);
            __m256 xx = _mm256_broadcast_ss(&inc);
            val = val + xx;
            _mm256_storeu_ps(dst, val);
            if(MAX_LOCAL > 8 && n>8){
                val = _mm256_loadu_ps(src+8);
                val = val + xx;
                _mm256_storeu_ps(dst+8, val);
            }
        }
    }

    struct polymer_bead_cache
    {
        // Spatial position
        vec3f_t x;
        // Which cell it is currently in
        uint32_t cell_index : 24;
        uint32_t cell_offset : 8;
    };

    struct Cell
    {
        unsigned cell_index;
        vec3i_t location;
        bool is_on_boundary;
        std::array<uint32_t,27> nhood_info; // Indexes of neighbouring cells. If is_on_boundary is true, then upper 6 bits say what wrapping to apply
        WorldState *world; // Convenience member for debug assertions

        // These are stored in structure of arrays form to make building full nhoods faster
        unsigned local_nhood_count;
        uint32_t ids[MAX_LOCAL];
        float positions[3][MAX_LOCAL];
        float velocities[3][MAX_LOCAL];        
        float forces[3][MAX_LOCAL];
        uint32_t bead_indices[MAX_LOCAL]; // Store the linear index rather than pointer to bead.

        uint8_t _pad_[64];

        // This set is shared between threads, and so can be contended
        std::atomic<unsigned> incoming_nhood_count;
        uint32_t incoming_ids[MAX_LOCAL];
        float incoming_positions[3][MAX_LOCAL];
        float incoming_velocities[3][MAX_LOCAL];        
        float incoming_forces[3][MAX_LOCAL];
        uint32_t incoming_bead_indices[MAX_LOCAL]; // Store the linear index rather than pointer to bead.

        tbb::concurrent_vector<Bead*> orphans;

        Cell()
            : local_nhood_count(0)
            , incoming_nhood_count(0)
        {}

        Cell(const Cell &x)
            : cell_index(x.cell_index)
            , location(x.location)
            , is_on_boundary(x.is_on_boundary)
            , nhood_info(x.nhood_info)
            , world(x.world)
            , local_nhood_count(0)
            , incoming_nhood_count(0)
        {
            if(x.local_nhood_count!=0){
                throw std::runtime_error("Only copy empty cells.");
            }
        }

        Cell &operator=(const Cell &o)
        {
            if(o.local_nhood_count!=0 || o.incoming_nhood_count!=0){
                throw std::runtime_error("Only copy empty cells.");
            }
            cell_index=o.cell_index;
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
        int add_to_local(uint32_t bead_index, uint32_t id, T position, T velocity, T force)
        {
            assert(world->beads[bead_index].get_hash_code().hash==id);

            if(local_nhood_count==MAX_LOCAL){
                return add_to_incoming(bead_index, id, position, velocity, force);
            }

            unsigned index=local_nhood_count++;
            bead_indices[index]=bead_index;
            ids[index]=id;
            for(int i=0; i<3; i++){
                positions[i][index] = position[i];
                velocities[i][index] = velocity[i];
                forces[i][index] = force[i];
            }

            return index;
        }

        template<class T>
        int add_to_incoming(uint32_t bead_index, uint32_t id, T position, T velocity, T force)
        {
            assert(world->beads[bead_index].get_hash_code().hash==id);

            unsigned index=incoming_nhood_count.fetch_add(1);
            if(index>=MAX_LOCAL){
                incoming_nhood_count.fetch_sub(1);
                return -1;
            }

            incoming_bead_indices[index]=bead_index;
            incoming_ids[index]=id;
            for(int i=0; i<3; i++){
                incoming_positions[i][index] = position[i];
                incoming_velocities[i][index] = velocity[i];
                incoming_forces[i][index] = force[i];
            }

            return index+MAX_LOCAL;
        }

        unsigned add_to_local(Bead *b)
        {
            return add_to_local(b->bead_id, b->get_hash_code().hash, b->x, b->v, b->f);
        }

        void move_incoming_to_local_and_sync_polymer_cache(std::vector<polymer_bead_cache> &cache)
        {
            unsigned nincoming=incoming_nhood_count.load(std::memory_order_relaxed);
            if(nincoming==0){
                return;
            }
         
            while(nincoming>0 && local_nhood_count<MAX_LOCAL){
                unsigned i = --nincoming;
                assert(  world->beads[incoming_bead_indices[i]].get_hash_code().hash == incoming_ids[i] );
                bead_indices[local_nhood_count]=incoming_bead_indices[i];
                ids[local_nhood_count]=incoming_ids[i];
                for(int d=0; d<3; d++){
                    positions[d][local_nhood_count]=incoming_positions[d][i];
                    velocities[d][local_nhood_count]=incoming_velocities[d][i];
                    forces[d][local_nhood_count]=incoming_forces[d][i];
                }

                if(!BeadHash{ids[local_nhood_count]}.is_monomer()){
                    auto &pc=cache[bead_indices[local_nhood_count]];
                    std::cerr<<" t="<<world->t<<", updating bead cache "<<bead_indices[i]<<"  (in incoming transfer to local)\n";
                    pc.cell_index=cell_index;
                    pc.cell_offset=local_nhood_count;
                    for(int d=0; d<3; d++){
                        pc.x[d]=positions[d][local_nhood_count];
                    }
                }

                ++local_nhood_count;
            }

            for(unsigned i=0; i<nincoming; i++){
                if(!BeadHash{incoming_ids[i]}.is_monomer()){
                    auto &pc=cache[incoming_bead_indices[i]];
                    std::cerr<<" t="<<world->t<<", updating bead cache "<<bead_indices[i]<<"  (in incoming transfer to non-local)\n";
                    pc.cell_index=cell_index;
                    pc.cell_offset=local_nhood_count+MAX_LOCAL;
                    for(int d=0; d<3; d++){
                        pc.x[d]=incoming_positions[d][local_nhood_count];
                    }
                }
            }

            incoming_nhood_count=nincoming;
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

        void sync_incoming_to_bead(WorldState *state, unsigned index)
        {
            assert(index < incoming_nhood_count);

            Bead &b=state->beads[incoming_bead_indices[index]];
            assert(b.bead_id==incoming_bead_indices[index]);
            assert(b.get_hash_code().hash==incoming_ids[index]);
            for(int i=0; i<3; i++){
                b.x[i] = incoming_positions[i][index];
                b.v[i] = incoming_velocities[i][index];
                b.f[i] = incoming_forces[i][index];
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

    void __attribute__((noinline)) import_beads_range(const range_t &r, tbb::concurrent_vector<const Bead*> &orphans)
    {
        for(unsigned i=r.begin(); i<r.end(); i++){
            const Bead &b=m_world->beads[i];
            vec3i_t cell_location=vec3_floor(b.x);
            unsigned cell_index=get_cell_index(cell_location);
            Cell &cell = m_cells[cell_index];

            // Pre-correct backwards by one step            
            auto velocity_half = b.v - b.f * (m_dt*0.5f);

            int cell_offset = cell.add_to_incoming( b.bead_id, b.get_hash_code().hash, b.x, velocity_half, b.f);
            if(cell_offset<0){
                orphans.push_back(&b);
            }

            if(!b.is_monomer){
                for(int d=0; d<3; d++){
                    m_polymer_positions[b.bead_id].x[d] = b.x[d];
                }
                m_polymer_positions[b.bead_id].cell_index=cell_index;
                m_polymer_positions[b.bead_id].cell_offset=cell_offset;
            }
        }
    }

    void __attribute__((noinline)) import_orphans_range(const range_t &r, tbb::concurrent_vector<const Bead*> &orphans)
    {
        for(unsigned i=r.begin(); i<r.end(); i++){
            const Bead &b=*orphans[i];
            vec3i_t cell_location=vec3_floor(b.x);
            unsigned cell_index=get_cell_index(cell_location);
            Cell &cell = m_cells[cell_index];

            // Pre-correct backwards by one step            
            auto velocity_half = b.v - b.f * (m_dt*0.5f);

            int cell_offset = cell.add_to_incoming( b.bead_id, b.get_hash_code().hash, b.x, velocity_half, b.f);
            if(cell_offset<0){
                throw std::runtime_error("MAX_LOCAL is too small");
            }

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
        tbb::concurrent_vector<const Bead*> orphans;
        orphans.reserve(m_world->beads.size());

        // First round we fill up incoming.

        ParallelFor<range_t>(range_t(0, m_world->beads.size(), 4096), [&](const range_t &r){
            import_beads_range(r, orphans);
        });

        ParallelFor<range_t>( range_t(0, m_cells.size(), 256), [&](const range_t &r){
            for(unsigned i=r.begin(); i<r.end(); i++){
                m_cells[i].move_incoming_to_local_and_sync_polymer_cache(m_polymer_positions);
            }
        });

        // Second round we deal out any orphans

        unsigned norphans=orphans.size();
        if(norphans){
            ParallelFor<range_t>(range_t(0, norphans, 1024), [&](const range_t &r){
                import_orphans_range(r, orphans);
            });

            ParallelFor<range_t>( range_t(0, m_cells.size(), 256), [&](const range_t &r){
                for(unsigned i=r.begin(); i<r.end(); i++){
                    m_cells[i].move_incoming_to_local_and_sync_polymer_cache(m_polymer_positions);
                }
            });
        }
    }

    void __attribute__((noinline)) export_beads_range(const range_t &r)
    {
        for(unsigned ci=r.begin(); ci<r.end(); ci++){
            Cell &c=m_cells[ci];
            for(unsigned i=0; i<c.local_nhood_count; i++){
                c.sync_local_to_bead(m_world, i);

                // Post-correct velocity forwards, as it hasn't been done
                auto &b=m_world->beads[c.bead_indices[i]];
                assert(!isnanf(b.f[0]));
                assert(!isnanf(b.v[0]));
                b.v = b.v + b.f * (m_dt*0.5f);
                assert(!isnanf(b.v[0]));
            }

            unsigned nincoming=c.incoming_nhood_count;
            for(unsigned i=0; i<nincoming; i++){
                c.sync_incoming_to_bead(m_world, i);

                // Post-correct velocity forwards, as it hasn't been done
                auto &b=m_world->beads[c.incoming_bead_indices[i]];
                assert(!isnanf(b.v[0]));
                b.v = b.v + b.f * (m_dt*0.5f);
                assert(!isnanf(b.v[0]));
            }

            // clear the cell out for the next round
            c.local_nhood_count=0;
            c.incoming_nhood_count=0;
        }
    }

    void export_beads()
    {
        static unsigned int EXPORT_BEADS_GRAIN = 256;

        ParallelFor(range_t(0, m_cells.size(), 256), [&](const range_t &r){
            export_beads_range(r);
        });
    }

    void validate_polymer_cache_range(const range_t &r)
    {
        for(unsigned i=r.begin(); i<r.end(); i++){
            Bead &b=m_world->beads[i];
            assert(b.bead_id==i);
            if(!b.is_monomer){
                auto pos=m_polymer_positions[i];
                const auto &cell=m_cells[pos.cell_index];
                if(pos.cell_offset < MAX_LOCAL){
                    assert( cell.bead_indices[pos.cell_offset] == i);
                    for(int d=0; d<3; d++){
                        if(cell.positions[d][pos.cell_offset] != pos.x[d] ){
                            std::cerr<<"  t="<<m_world->t<<", id="<<i<<", true={"<<cell.positions[0][pos.cell_offset]<<","<<cell.positions[1][pos.cell_offset]<<","<<cell.positions[2][pos.cell_offset]<<"}, got={"<<pos.x[0]<<","<<pos.x[1]<<","<<pos.x[2]<<"}\n";
                        }
                        assert( cell.positions[d][pos.cell_offset] == pos.x[d] );
                    }
                }else{
                    assert( cell.bead_indices[pos.cell_offset-MAX_LOCAL] == i);
                    for(int d=0; d<3; d++){
                        assert( cell.positions[d][pos.cell_offset-MAX_LOCAL] == pos.x[d] );
                    }
                }
            }
        }
    }

    void validate_polymer_cache()
    {
        ParallelFor(range_t(0, m_world->beads.size(), 8192), [&](const range_t &r){
            validate_polymer_cache_range(r);
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

        unsigned total=cell.local_nhood_count + cell.incoming_nhood_count.load();
        for(unsigned local_index=0; local_index<cell.local_nhood_count; local_index++){
            float position[3], velocity[3], force[3];
            uint32_t id;

            if( local_index < MAX_LOCAL ){
                for(unsigned d=0; d<3; d++){
                    position[d]=cell.positions[d][local_index];
                    velocity[d]=cell.velocities[d][local_index];
                    force[d]=0.0f;
                }
                id=cell.ids[local_index];
            }else{
                for(unsigned d=0; d<3; d++){
                    position[d]=cell.incoming_positions[d][local_index-MAX_LOCAL];
                    velocity[d]=cell.incoming_velocities[d][local_index-MAX_LOCAL];
                    force[d]=0.0f;
                }
                id=cell.incoming_ids[local_index-MAX_LOCAL];
            }


            unsigned bead_type=BeadHash{id}.get_bead_type();

            // We have to do this out of the main loop, as it requires things
            // that gcc can't inferr (shuffle)
            const float *conStrengthRow = conStrengthMatrix+bead_type*MAX_BEAD_TYPES;
            const float *sqrtDissStrengthRow = sqrtDissStrengthMatrix+bead_type*MAX_BEAD_TYPES;
            if( (MAX_LOCAL==8 || MAX_LOCAL==16) && MAX_BEAD_TYPES <= 8){
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
            }

            if(local_index<MAX_LOCAL){
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    cell.forces[d][local_index] += force[d];
                }
            }else{
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    cell.incoming_forces[d][local_index-MAX_LOCAL] += force[d];
                }
            }
        }
    }

    void __attribute__((noinline)) move_beads_range(const range_t &r)
    {
        float dt=m_dt;  

        for(unsigned index=r.begin(); index<r.end(); index++){
            auto &cell=m_cells[index];              // Exclusive access
            auto &cell_next=m_cells_next[index];    // Exclusive access to local (others can access incoming)

            /*if(index+1<r.end()){
                _mm_prefetch( m_cells[index+1].bead_indices, _MM_HINT_T0);
                _mm_prefetch( m_cells[index+1].forces, _MM_HINT_T0 );
                _mm_prefetch( m_cells[index+1].positions, _MM_HINT_T0 );
                _mm_prefetch( m_cells[index+1].velocities, _MM_HINT_T0 );
            }*/

            // We use this to handle both normal and incoming with the same code.
            // It is rare we come back.
            bool take_slow_path=false;
        repeat_movement:

            unsigned n=cell.local_nhood_count;
            declare(n <= MAX_LOCAL);

            bool has_moved[MAX_LOCAL];
            bool is_polymer[MAX_LOCAL];
            int32_t dest_cell_index[3][MAX_LOCAL];
            bool any_moved=0;
            bool any_polymer=false;

            float position[3][MAX_LOCAL];
            float velocity[3][MAX_LOCAL];
            float force[3][MAX_LOCAL];
            for(unsigned i=0; i<n; i++){
                assert(  m_world->beads[cell.bead_indices[i]].get_hash_code().hash == cell.ids[i] );
                
                uint32_t bead_index;
                uint32_t id;

                id=cell.ids[i];
                bead_index=cell.bead_indices[i];

                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    position[d][i]=cell.positions[d][i];
                    velocity[d][i]=cell.velocities[d][i];
                    force[d][i]=cell.forces[d][i];
                }

                // Mom step to get from v(t-dt/2) to v(t) (needs to be pre-corrected before first step)
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    velocity[d][i] += force[d][i] * (dt*0.5f);
                }

                // x(t+dt) = x(t) + v(t)*dt + f(t)*dt*dt/2
                // v(t+dt/2) = v(t) + dt*f(t)/2
                
                has_moved[i]=0;
                #pragma GCC unroll 3
                for(int d=0; d<3; d++){
                    position[d][i] += velocity[d][i]*dt + force[d][i] * dt*(dt*0.5f);
                    velocity[d][i] += force[d][i] * (dt*0.5f);
                    force[d][i] = 0;
                    dest_cell_index[d][i] = floorf(position[d][i]);
                    has_moved[i] |= dest_cell_index[d][i] != cell.location[d]; 
                }
                any_moved |= has_moved[i];
                any_polymer |= !BeadHash{id}.is_monomer();

                assert(  m_world->beads[cell.bead_indices[i]].get_hash_code().hash == cell.ids[i] );
            }

            // We want to avoid the fast path if we have to repeat with incoming,
            // as the relative indices will be messed up through the fast-path copying.
            if(!any_moved && !take_slow_path){
                // Hopefully... the common + fast path

                cell_next.local_nhood_count = n;
                copy_max_local_words(cell_next.ids, cell.ids, n);
                copy_max_local_words(cell_next.bead_indices, cell.bead_indices, n);
                for(unsigned d=0; d<3; d++){
                    copy_max_local_words(cell_next.positions[d], position[d], n);
                    copy_max_local_words(cell_next.velocities[d], velocity[d], n);
                    copy_max_local_words(cell_next.forces[d], force[d], n);
                }

                if(any_polymer){
                    // Comparitively rare
                    for(unsigned i=0; i<n; i++){
                        if(!BeadHash{cell.ids[i]}.is_monomer()){
                            uint32_t bead_index=cell.bead_indices[i];
                            std::cerr<<" t="<<m_world->t<<", updating bead cache "<<i<<"  (in fast path)\n";
                            for(unsigned d=0; d<3; d++){
                                m_polymer_positions[bead_index].x[d] = position[d][i];
                            }
                            m_polymer_positions[bead_index].cell_index=cell.cell_index;
                            m_polymer_positions[bead_index].cell_offset=i;
                        }
                    }
                }
            }else{
                // The slow path. Hopefully beads leaving is relatively rare (e.g. 1-10% / step)
                for(unsigned i=0; i<n; i++){
                    // Sigh.
                    float lpos[3]={position[0][i], position[1][i], position[2][i] };
                    float lvel[3]={velocity[0][i], velocity[1][i], velocity[2][i] };
                    float lfor[3]={force[0][i], force[1][i], force[2][i] };
                    if(has_moved[i]){
                        for(unsigned d=0; d<3; d++){
                            if(dest_cell_index[d][i] < 0){
                                assert(lpos[d] >= -1.0f);
                                lpos[d] += m_dimsf[d];
                                dest_cell_index[d][i] += m_dims[d];
                                if(lpos[d] == m_dimsf[d]){ // Stupid edge case when position is _just_ under 0
                                    assert(dest_cell_index[d][i]==m_dims[i]-1);
                                    lpos[d] = 0;
                                    dest_cell_index[d][i] = 0;
                                }
                            }else if(dest_cell_index[d][i] >= m_dims[d]){
                                lpos[d] -= m_dimsf[d];
                                dest_cell_index[d][i] -= m_dims[d];
                            }
                        }

                        assert(  m_world->beads[cell.bead_indices[i]].get_hash_code().hash == cell.ids[i] );
                        unsigned dst_index=get_cell_index(vec3i_t{dest_cell_index[0][i], dest_cell_index[1][i], dest_cell_index[2][i]});
                        if(0 > m_cells_next.at(dst_index).add_to_incoming(cell.bead_indices[i], cell.ids[i], lpos, lvel, lfor)){
                            throw std::runtime_error("MAX_LOCAL is too small");
                        }

                        // any polymer will be synced while moving from incoming to local
                    }else{
                        int cell_offset = cell_next.add_to_local(cell.bead_indices[i], cell.ids[i], lpos, lvel, lfor);

                        if(is_polymer[i]){
                            std::cerr<<" t="<<m_world->t<<", updating bead cache "<<i<<"  (in slow local path)\n";
                            uint32_t bead_index=cell.bead_indices[i];
                            for(int d=0; d<3; d++){
                                m_polymer_positions[bead_index].x[d] = position[d][i];
                            }
                            m_polymer_positions[bead_index].cell_index=cell.cell_index;
                            m_polymer_positions[bead_index].cell_offset=cell_offset;
                        }
                    }
                }
            }

            cell.local_nhood_count=0;

            unsigned nincoming=cell.incoming_nhood_count;
            if(nincoming){
                cell.local_nhood_count=nincoming;
                copy_max_local_words(cell.ids, cell.incoming_ids, nincoming);
                copy_max_local_words(cell.bead_indices, cell.incoming_bead_indices, nincoming);
                for(int d=0; d<3; d++){
                    copy_max_local_words(cell.positions[d], cell.incoming_positions[d], nincoming);
                    copy_max_local_words(cell.velocities[d], cell.incoming_velocities[d], nincoming);
                    copy_max_local_words(cell.forces[d], cell.incoming_forces[d], nincoming);
                }

                take_slow_path=true;
                cell.incoming_nhood_count=0;
                goto repeat_movement;
            }
        }
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
                    (head_pos.cell_offset>=MAX_LOCAL ? head_cell.forces[d][head_pos.cell_offset-MAX_LOCAL] : head_cell.forces[d][head_pos.cell_offset] ) += force_head[d];
                    (tail_pos.cell_offset>=MAX_LOCAL ? tail_cell.forces[d][tail_pos.cell_offset-MAX_LOCAL] : tail_cell.forces[d][tail_pos.cell_offset] ) -= force_head[d];
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
                for(int d=0; d<3; d++){
                    if(head_pos.cell_offset<MAX_LOCAL){
                        assert(m_cells[head_pos.cell_index].positions[d][head_pos.cell_offset] == head_pos.x[d] );
                    }else{
                        assert(m_cells[head_pos.cell_index].incoming_positions[d][head_pos.cell_offset-MAX_LOCAL] == head_pos.x[d] );
                    }
                }
                auto centre_pos=m_polymer_positions[centre_index];
                auto tail_pos=m_polymer_positions[tail_index];

                vec3f_t dx01 = centre_pos.x - head_pos.x ;
                vec3f_t dx12 = tail_pos.x - centre_pos.x;
                for(int d=0; d<3; d++){
                    if(dx01[d] > 2.5){
                        dx01[d] -= m_dimsf[d];
                    }else if(dx01[d] < -2.5){
                        dx01[d] += m_dimsf[d];
                    }
                    if(dx12[d] > 2.5){
                        dx12[d] -= m_dimsf[d];
                    }else if(dx12[d] < -2.5){
                        dx12[d] += m_dimsf[d];
                    }
                }

                assert(dx01.l2_norm()!=0);
                assert(dx12.l2_norm()!=0);
                assert(dx01.l2_norm()<1);
                assert(dx12.l2_norm()<1);

                // TODO : Calculating this every time is stupid
                float sin_theta, cos_theta;
                sincosf(bp.theta0, &sin_theta, &cos_theta);

                vec3f_t head_force, centre_force, tail_force;
                dpd_maths_core::calc_angle_force(
                    float(bp.kappa), cos_theta, sin_theta,
                    dx01, dx01.l2_norm(), dx12, dx12.l2_norm(),
                    head_force, centre_force, tail_force
                );
                /*if(isnanf(head_force[0])){
                    std::cerr<<"dx01="<<dx01<<", dx12"<<dx12<<"\n";
                }*/

                //std::cerr<<"("<<bp.bond_offset_head<<","<<bp.bond_offset_tail<<"), dx01="<<dx01<<", dx12"<<dx12<<"\n";

                assert(!isnanf(head_force[0]));
                assert(!isnanf(centre_force[0]));
                assert(!isnanf(tail_force[0]));

                Cell &head_cell=m_cells[head_pos.cell_index];
                Cell &centre_cell=m_cells[centre_pos.cell_index];
                Cell &tail_cell=m_cells[tail_pos.cell_index];

                for(int d=0; d<3; d++){
                    (head_pos.cell_offset>=MAX_LOCAL ? head_cell.forces[d][head_pos.cell_offset-MAX_LOCAL] : head_cell.forces[d][head_pos.cell_offset] ) += head_force[d];
                    (centre_pos.cell_offset>=MAX_LOCAL ? centre_cell.forces[d][centre_pos.cell_offset-MAX_LOCAL] : centre_cell.forces[d][centre_pos.cell_offset] ) += centre_force[d];
                    (tail_pos.cell_offset>=MAX_LOCAL ? tail_cell.forces[d][tail_pos.cell_offset-MAX_LOCAL] : tail_cell.forces[d][tail_pos.cell_offset] ) += tail_force[d];
                }

                /*if(ForceLogging::logger()){
                    auto h1=m_world->beads[head_index].get_hash_code();
                    auto h2=m_world->beads[centre_index].get_hash_code();
                    auto h3=m_world->beads[tail_index].get_hash_code();
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_xhd", head_pos.x);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_xmd", centre_pos.x);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_xtl", tail_pos.x);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_dx01", dx01);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_dx12", dx12);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_head", head_force);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_mid", centre_force);
                    ForceLogging::logger()->LogBeadTripleProperty(h1, h2, h3, "f_next_angle_tail", tail_force);
                }*/
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

        validate_polymer_cache();

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

        /*for(auto &c : m_cells){
            for(unsigned i=0; i<c.local_nhood_count; i++){
                assert(  m_world->beads[c.bead_indices[i]].get_hash_code().hash == c.ids[i] );
            }
        }*/

        ParallelFor<range_t>( range_t(0, m_cells.size(), 256), [&](const range_t &r){
            move_beads_range(r);
        });

        std::swap(m_cells, m_cells_next);

        /*for(auto &c : m_cells){
            for(unsigned i=0; i<c.local_nhood_count; i++){
                assert(  m_world->beads[c.bead_indices[i]].get_hash_code().hash == c.ids[i] );
            }
            for(unsigned i=0; i<c.incoming_nhood_count; i++){
                assert(  m_world->beads[c.incoming_bead_indices[i]].get_hash_code().hash == c.incoming_ids[i] );
            }
        }*/

        ParallelFor<range_t>( range_t(0, m_cells.size(), 1024), [&](const range_t &r){
            for(unsigned i=r.begin(); i<r.end(); i++){
                m_cells[i].move_incoming_to_local_and_sync_polymer_cache(m_polymer_positions);
            }
        });

        validate_polymer_cache();

        /*for(auto &c : m_cells){
            for(unsigned i=0; i<c.local_nhood_count; i++){
                assert(  m_world->beads[c.bead_indices[i]].get_hash_code().hash == c.ids[i] );
            }
        }*/

        ParallelFor<range_t>( range_t(0, m_world->polymers.size(), 64), [&](const range_t &r){
            update_polymer_range(r);
        });

        // We want a reasonable grain size due to the allocation for working set and
        // the caching opportunities between adjacent cells. Maybe this should be a 3d
        // range?
        ParallelFor<range_t>( range_t(0, m_cells.size(), 256), [&](const range_t &r){
            update_cell_force_range(r);
        }); 

        /*for(auto &c : m_cells){
            for(unsigned i=0; i<c.local_nhood_count; i++){
                assert(  m_world->beads[c.bead_indices[i]].get_hash_code().hash == c.ids[i] );
            }
        }*/
    }
public:


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

        ParallelFor(range_t(0, ncells, 64), [&](const range_t &r){
            for(unsigned i=r.begin(); i<r.end(); i++){
                vec3i_t location=get_cell_pos(i);
                Cell &c=m_cells[i];
                c.cell_index=i;
                c.location=location;
                c.is_on_boundary=false;
                for(int d=0; d<3; d++){
                    c.is_on_boundary |= (location[d]==0) || (location[d]==m_dims[d]-1);
                }
                c.world=m_world;
                c.local_nhood_count=0;
                c.incoming_nhood_count=0;

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
            }
        });
    }

    bool CanSupportHookeanBonds() const
    { return true; }

    bool CanSupportAngleBonds() const
    { return true; } // Let the specific error messagecome through.

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

        return DPDEngine::CanSupport(world);
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
