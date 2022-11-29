#include <cstdint>
#include <cmath>
#include <cassert>

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace dpd_maths_highway
{
    namespace HWY_NAMESPACE
    {
        constexpr int MAX_BEAD_TYPES = 8;

        namespace hn = hwy::HWY_NAMESPACE;

        using vecf_tag = hn::ScalableTag<float>;

        constexpr int K = MaxLanes(vecf_tag());

        using vecu_tag = hn::ScalableTag<uint32_t>;
        using vecs_tag = hn::ScalableTag<int32_t>;
        static_assert( MaxLanes(vecf_tag())==MaxLanes(vecu_tag()) );

        struct packed_bead
        {
            uint32_t hash;
            float x[3];
            float v[3];
            float f[3];
        };

        /*
        Structure optimised for minimum size vectorised access and force accumulation
        This assumes that unaligned loads and stores are fast, and that total cache
        traffic is the more problematic thing to deal with.
        For example, if we have 64-byte / 16 word cache line:
        - soa_packed_beads, stride=3: (1+10*3) = 31 words -> 2 cache lines
        - soa_packed_beads, stride=4: (1+10*4) = 41 words -> 3 cache lines
        - soa avx2 : 1+10*8 = 81 words -> 6 cache lines
        - soa avx512 : 1+10*16 = 161 words -> 11 cache lines

        It's a pain to deal with when adding/removing beads, but given most accesses
        are read this is ok. The stride parameter means we can have n<=stride, so
        in many cases we could keep stride at some parameter and have small gaps.
        Occasionally it can be re-packed to minimise the number of cache blocks.
        */
        struct soa_packed_beads
        {
            static constexpr int MAX_BEADS = 16;

            /*
            Each bead has:
            - 32-bit hash
            - 3-vec float position
            - 3-vec float velocity
            - 3-vec float force
            */
           uint32_t pad_lo[14]; // Allows us to perform unaligned reads before the true data
           uint32_t stride;     // Length of vectors
           uint32_t n;          // Number of valid beads
           uint32_t data[10*MAX_BEADS];
           uint32_t pad_hi[16]; // Allows us to perform unaligned reads and writes to the force vector

           void clear(uint32_t nstride=0)
           {
            assert(nstride <= MAX_BEADS);
            stride=nstride;
            n=0;
           }

           void set_stride(uint32_t nstride)
           {
            assert(nstride >= n);
            assert(nstride <= MAX_BEADS);
            if(n==0){
                stride=nstride;
            }else if(nstride>stride){
                for(int i=9; i>0; i--){
                    std::copy(data+stride*i, data+stride*i+n, data+nstride*i);
                }
            }else if(nstride<stride){
                for(int i=1; i<=9; i++){
                    std::copy_backward(data+stride*i, data+stride*i+n, data+nstride*i);
                }
            }
           }

           packed_bead erase(unsigned index)
           {
            float *fdata=(float*)data;

            assert(index < n);
            packed_bead res;
            res.hash=data[index];
            for(int d=0; d<3; d++){
                res.x[d]=fdata[ (1+d)*stride ];
                res.v[d]=fdata[ (4+d)*stride ];
                res.f[d]=fdata[ (7+d)*stride ];
            }
            n -= 1;
            if(index != n){
                data[index] = data[n];
                for(int d=0; d<3; d++){
                    fdata[ (1+d)*index ] = fdata[ (1+d)*n ];
                    fdata[ (4+d)*index ] = fdata[ (4+d)*n ];
                    fdata[ (7+d)*index ] = fdata[ (7+d)*n ];
                }                
            }
            return res;
           }

           void insert(const packed_bead &b)
           {
            assert(n < MAX_BEADS);
            if(n == stride){
                set_stride(stride+1);
            }
            assert(n<stride);
            float *fdata=(float*)data;
            unsigned index=n;

            data[index] = b.hash;
            for(int d=0; d<3; d++){
                fdata[ (1+d)*index ] = b.x[d];
                fdata[ (4+d)*index ] = b.v[d];
                fdata[ (7+d)*index ] = b.f[d];
            }
           }

           const uint32_t * HWY_RESTRICT get_hash_vec() const
           { return data+0; } 

           const float * HWY_RESTRICT get_x_vec(int d) const
           { return (const float*)data+(1+d)*stride; }

           const float * HWY_RESTRICT get_v_vec(int d) const
           { return (const float*)data+(4+d)*stride; } 

           const float * HWY_RESTRICT get_f_vec(int d) const
           { return (const float*)data+(7+d)*stride; } 

           float * HWY_RESTRICT get_f_vec(int d)
           { return (float*)data+(7+d)*stride; } 

           void newton(float origin[3], float half_dt, std::vector<uint8_t> &outgoing)
           {
            outgoing.clear();

            int offset=0;
            while(off < n){
                hn::Vec<vecf_tag> x[3], v[3], f[3];
                hn::Mask<vecf_tag> moved;
                for(int d=0; d<3; d++){
                    x[d]=hn::LoadU(vecf_tag(), get_x_vec(d)+offset);
                    v[d]=hn::loadU(vecf_tag(), get_v_vec(d)+offset);
                    f[d]=hn::loadU(vecf_tag(), get_f_vec(d)+offset);
                
                    // update mom
                    Vec<vecf_tag> half_dv = hn::Set(vecf_tag(), half_dt) * f[d];
                    v[d] += half_dv;

                    // Update x
                    // auto x = b.x + b.v*dt + b.f*half(dt*dt);
                    x[d] += hn::Set(vecf_tag(), dt) * ( v[d] + f[d]*hn::Set(vecf_tag(), half_dt));
                    v[d] += half_dv;

                    auto dim_moved = x[d] < hn::Set(vecf_tag(), origin[d]) || x[d] >= hn::Set(vecf_tag(), origin[d]+1.0f);
                    if(d==0){
                        moved=dim_moved;
                    }else{
                        moved |= dim_moved;
                    }
                }

                // Must be done in sequence, as they may overlap
                hn::StoreU( x[0], vecf_tag(), get_x_vec(0)+offset );
                hn::StoreU( x[1], vecf_tag(), get_x_vec(1)+offset );
                hn::StoreU( x[2], vecf_tag(), get_x_vec(2)+offset );
                hn::StoreU( v[0], vecf_tag(), get_v_vec(0)+offset );
                hn::StoreU( v[1], vecf_tag(), get_v_vec(1)+offset );
                hn::StoreU( v[2], vecf_tag(), get_v_vec(2)+offset );
                hn::StoreU( f[0], vecf_tag(), get_f_vec(0)+offset );
                hn::StoreU( f[1], vecf_tag(), get_f_vec(1)+offset );
                hn::StoreU( f[2], vecf_tag(), get_f_vec(2)+offset );

                if(!hn::AllFalse(moved)){
                    int f=hn::FindFirstTrue(moved);
                    while(f!=-1){
                        outgoing.push_back(f);
                        
                        moved = hn::And( moved, hn::Not(hn::FirstN(vecf_tag(), f+1)) );
                        f=hn::FindFirstTrue(moved);
                    }
                }

                offset += K;
            }
           }
        };

        struct vector_packed_beads
        {
            static constexpr int MAX_BEADS = 16 * soa_packed_beads::MAX_BEADS;

            unsigned n;

            uint32_t hash[MAX_BEADS];
            float x[3][MAX_BEADS];
            float v[3][MAX_BEADS];
            float f[3][MAX_BEADS];

            void gather_beads(
                unsigned ncells,
                const soa_packed_beads *cells
            ){
                unsigned offset = 0;
                
                n=0;
                for(unsigned i=0; i<ncells; i++){
                    const soa_packed_beads &cell = cells[i];
                    int off=0;
                    while(off < cell.n){
                        hn::StoreU( hn::LoadU( vecu_tag(), cell.get_hash_vec()+off), vecu_tag(), hash+n );
                        for(int d=0; d<3; d++){
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_x_vec(d)+off), vecf_tag(), x[d]+n );
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_v_vec(d)+off), vecf_tag(), v[d]+n );
                        }
                        
                        n += std::min<unsigned>(cell.n-off, K );
                        off += K;
                    }

                    n += cell.n;
                }

                // Zero out all forces.
                // Hopefully all cache hits.
                for(int d=0; d<3; d++){
                    memset(f[d], 0, ((n+K-1) & (1-K))*4); 
                }
            }

            /*
            Adds the forces in the cells.
            Assumes we have exclusive access to destinations.
            */
            void scatter_forces(
                unsigned ncells,
                soa_packed_beads *cells
            ){
                unsigned offset=0;
                for(unsigned i=0; i<ncells; i++){
                    soa_packed_beads &cell=cells[i];

                    // Have to do in strict dimension order, as earlier dims
                    // clobber later dims. Better to do gather of forces into
                    // contiguous array, then add and store it packed
                    float local[3*soa_packed_beads::MAX_BEADS];
                    for(int d=0; d<3; d++){
                        std::copy( f[d]+offset, f[d]+offset+cell.n, local+cell.n*d );
                    }
                    wedge=local;
                    float *dst=cell.get_f_vec(0);
                    for(int v=0; v<cell.n; v+= K){
                        auto f_acc=hn::LoadU( vecf_tag(), dst );
                        auto f_new=hn::LoadU( vecf_tag(), wedge );
                        hn::StoreU( f_acc+f_new, vecf_tag(), dst );
                        dst += K;
                        wedge += K;
                    }

                    offset += cell.n;
                }
            }
        };


        /* Takes
        A set of packed beads and performs all interations.
        "Home" beads are in [0,nHome)
        All beads are in [0,nTotal)
        Broadly this does:
            for hi in [0,nHome):
                for oi in [hi+1,nTotal):
                    f=interact(beads[hi], beads[oi])
                    beads[hi].f += f
                    beads[lo].f -= f
        */
         void filter_and_calc_force(
            float rng_scale_s32_to_u, // Converts a value in [-2^31,2^31) to correctly scaled float
            uint32_t t_hash,

            const float conStrength[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
            const float sqrtDissStrength[MAX_BEAD_TYPES*MAX_BEAD_TYPES],

            unsigned nTotal,
            unsigned nHome,  // The first nHome are active interactors

            const uint32_t * HWY_RESTRICT hash, // top 4 bits are bead type index
            const float * HWY_RESTRICT x[3],
            const float * HWY_RESTRICT v[3],
            float * HWY_RESTRICT f[3]
        ) {
            unsigned nTotalVectors = (nTotal + Lanes(vecf_tag()) - 1) / Lanes(vecf_tag());
            unsigned nTotalPadded = nTotalVectors * Lanes(vecf_tag());

            for(unsigned hi=0; hi<nHome; hi++){
                float hx[3], hv[3];
                uint32_t hhash;
                for(unsigned d=0; d<3; d++){
                    hx[d] = x[d][hi];
                    hv[d] = v[d][hi];
                }
                hhash=hash[hi];
                hn::Vec<vecf_tag> conStrengthTable, sqrtDissStrengthTable;
                const float * HWY_RESTRICT conStrengthRow = conStrength+(hhash>>(32-4))*MAX_BEAD_TYPES;
                const float * HWY_RESTRICT sqrtDissStrengthRow=sqrtDissStrength+(hhash>>(32-4))*MAX_BEAD_TYPES;
                if constexpr( MAX_BEAD_TYPES <= MaxLanes(vecf_tag()) ){
                    conStrengthTable = hn::Load(vecf_tag(), conStrengthRow);
                    sqrtDissStrengthTable=hn::Load(vecf_tag(), sqrtDissStrengthRow);
                }

                hn::Vec<vecf_tag> f_home_acc[3]={ hn::Zero(vecu_tag()), hn::Zero(vecu_tag()), hn::Zero(vecu_tag()) };
                for(unsigned oi; oi<nTotalPadded; oi+=hn::Lanes(vecf_tag())){
                    hn::Vec<vecf_tag> ox[3], dx[3];
                    for(int d=0; d<3; d++){
                        ox[d]=hn::Load(vecf_tag(), x[d]+oi);
                        dx[d]=hn::Sub(hn::Set(vecf_tag(), hx[d]) , ox[d]);
                    }
                    auto dr2=dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                    if(0 && hn::AllTrue(dr2 >= 1)){
                        // TODO: When does this make sense? Probably never, even if K==4
                        // Maybe if the corners are placed together?
                        continue; 
                    }

                    auto ohash=hn::Load(vecu_tag(), hash+oi);
                    auto obead_type= hn::ShiftRight<32-4>(ohash);

                    hn::Vec<vecf_tag> ov[3], dv[3];
                    for(int d=0; d<3; d++){
                        ov[d]=hn::Load(vecf_tag(), v[d]+oi);
                        dv[d]=hn::Sub(hn::Set(vecf_tag(), hv[d]), ov[d]);
                    }
                    auto dr=hn::Sqrt(dr2);
                    const auto one=hn::Set(vecf_tag(), 1.0f);
                    auto wr= one - dr;
                    auto inv_dr= one/dr;

                    hn::Vec<vecf_tag> conStrengthVal, sqrtDissStrengthVal;
                    if constexpr( MAX_BEAD_TYPES <= MaxLanes(vecf_tag()) ){
                        auto index = hn::IndicesFromVec(vecf_tag(), obead_type);
                        conStrengthVal = hn::TableLookupLanes(conStrengthTable, index);
                        sqrtDissStrengthVal = hn::TableLookupLanes(sqrtDissStrengthTable, index);
                    }else{
                        conStrengthVal = hn::GatherIndex(vecf_tag(), conStrengthRow, hn::BitCast(vecs_tag(), obead_type));
                        sqrtDissStrengthVal = hn::GatherIndex(vecf_tag(), sqrtDissStrengthRow, hn::BitCast(vecs_tag(), obead_type));
                    }
                    
                    auto conForce = conStrengthVal*wr;

                    // dx is unnormalised
                    auto rdotv_mul_dr = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

                    auto sqrt_gammap = sqrtDissStrengthVal*wr;
                    auto negDissForce = sqrt_gammap*sqrt_gammap*rdotv_mul_dr*inv_dr;

                    auto ubits = hn::Set(vecu_tag(), hhash + t_hash) ^ (ohash + hn::Set(vecu_tag(), t_hash)); // TODO: fix
                    auto ubits_s = hn::BitCast(vecs_tag(), ubits);
                    auto u = hn::ConvertTo(vecf_tag(), ubits_s) * hn::Set(vecf_tag(), rng_scale_s32_to_u ); // TODO: Check range
                
                    auto randForce = sqrt_gammap * u;

                    auto scaled_force = (conForce - negDissForce + randForce ) * inv_dr;
                    auto valid_r2_lt_1 = dr > hn::Set(vecf_tag(), 0.00001f);
                    auto valid_skip_done = hn::Set(vecu_tag(), hi) < hn::Iota(uvec_tag(), oi);
                    scaled_force = hn::IfThenElseZero( valid_r2_lt_1 && valid_skip_done, scaled_force);

                    for(int d=0; d<3; d++){
                        auto f_home = scaled_force*dx[d];
                        f_home_acc[d] += f_home;
                        
                        auto ff_o = hn::Load(vecf_tag(), f[d]+oi);
                        // dx is unnormalised, but scale_force incorporates inv_dr
                        auto ff_n=ff - f_home;
                        float * HWY_RESTRICT ff_dst = f[d]+oi;
                        hn::Store(ff_n, vecf_tag(), ff_dst);
                    }
                }

                for(int d=0; d<3; d++){
                    f[d][hi] += hn::SumOfLanes(f_home_acc[d]);
                }
            }
        };


        void filter_and_calc_force(
            float rng_scale_s32_to_u,
            uint32_t t_hash,

            const float conStrength[MAX_BEAD_TYPES*MAX_BEAD_TYPES],
            const float sqrtDissStrength[MAX_BEAD_TYPES*MAX_BEAD_TYPES],

            unsigned ncells,
            soa_packed_beads *cells,

            vector_packed_beads *working
        ){
            working->gather_beads(ncells, cells);

            filter_and_calc_force(
                rng_scale_s32_to_u,
                t_hash,
                
                conStrength,
                sqrtDissStrength,

                working->n,
                cells[0].n,
                
                working->hash,
                working->x,
                working->v,
                working->f
            );


            working->scatter_forces(ncells, cells);

        }
    }

};
HWY_AFTER_NAMESPACE();
