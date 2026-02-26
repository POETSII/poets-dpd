#ifndef dpd_maths_highway_hpp
#define dpd_maths_highway_hpp

#include <cstdint>
#include <cmath>
#include <cassert>

#include "dpd/maths/dpd_maths_core.hpp"
#include "dpd/core/logging.hpp"

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace dpd_maths_highway
{
    struct packed_bead
    {
        uint32_t hash;
        float x[3];
        float v[3];
        float f[3];
    };

    struct ExecConfig
    {
        static constexpr int MAX_BEAD_TYPES = 8;
        static constexpr bool SHARED_SQRT_DISS = false;
        static constexpr bool USE_HASH = true;
    };


    namespace HWY_NAMESPACE
    {

        namespace hn = hwy::HWY_NAMESPACE;

        using vecf_tag = hn::ScalableTag<float>;

        constexpr int K = MaxLanes(vecf_tag());

        using vecu_tag = hn::ScalableTag<uint32_t>;
        using vecs_tag = hn::ScalableTag<int32_t>;
        static_assert( MaxLanes(vecf_tag())==MaxLanes(vecu_tag()) );

        // We have 26 neighbours plus self
        // Using hemi sharing we have 13 + 1
        static constexpr int NHOOD_SIZE = 14;

        struct LaneRotator
        {
            template<class D>
            using rot_type = decltype( hn::IndicesFromVec(D(), hn::Zero(D()) ) );

            using rot_token = rot_type<vecu_tag>;

            rot_type<vecu_tag> table[K];

            template<class D>
            rot_type<D> get_rot_impl(unsigned distance)
            {
                assert(distance<K);
                static const auto MASK=hn::Set(D(), K-1);
                static const auto IOTA_PLUS_K=hn::Set(D(), K);
                auto srcs= hn::And( hn::Sub(IOTA_PLUS_K,hn::Set(D(), distance)) , MASK);
                return hn::IndicesFromVec(D(), srcs);
            }

            LaneRotator()
            {
                for(unsigned i=0; i<K; i++){
                    table[i] = get_rot_impl<vecu_tag>(i);
                }
            }

            rot_token get_rot(unsigned distance)
            {
                return table[distance];
            }

            template<class D, class T>
            T apply_rot(D d, T x, rot_token rot)
            {
                using rot_type_tgt = decltype( hn::IndicesFromVec(d, hn::Zero(d) ) );

                // There is no safe way of rebinding in highway I'm away of
                static_assert(sizeof(rot_token)==sizeof(rot_type_tgt));

                auto trot =reinterpret_cast<rot_type_tgt>(rot);

                return hn::TableLookupLanes( x, trot );
            }
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
           hn::Vec<vecu_tag> _pad_align_;
           uint32_t _align_lo_[128-sizeof(_pad_align_)-2*4];
           uint32_t stride;     // Length of vectors
           uint32_t n;          // Number of valid beads
           uint32_t data[10*MAX_BEADS];
           uint32_t pad_hi[16]; // Allows us to perform unaligned reads and writes to the force vector

            soa_packed_beads()
                : stride(0)
                , n(0)
            {}

           void clear(uint32_t nstride=0)
           {
            assert(nstride <= MAX_BEADS);
            stride=nstride;
            n=0;
           }

           void set_stride(uint32_t nstride)
           {    
            #if 0
            std::vector<packed_bead> prev;
            for(unsigned i=0; i<n; i++){
                prev.push_back(get_bead(i));
            }
            n=0;
            stride=nstride;
            for(unsigned i=0; i<prev.size(); i++){
                insert(prev[i]);
                auto g=get_bead(i);
                assert(!memcmp(&g, &prev[i], sizeof(packed_bead)));
            }

            #else


            /*#ifndef NDEBUG
            packed_bead gg;
            if(n>0){
                gg=get_bead(0);
            }
            std::vector<packed_bead> prev;
            for(unsigned i=0; i<n; i++){
                prev.push_back(get_bead(i));
            }
            #endif*/

            assert(nstride >= n);
            assert(nstride <= MAX_BEADS);
            if(n==0){
                // do nothing
            }else if(nstride>stride){
                for(int i=9; i>0; i--){
                    assert(data+stride*i < data+nstride*i);
                    std::copy_backward(data+stride*i, data+stride*i+n, data+nstride*i+n);
                }
            }else if(nstride<stride){
                for(int i=1; i<=9; i++){
                    assert(data+stride*i > data+nstride*i);
                    std::copy(data+stride*i, data+stride*i+n, data+nstride*i);
                }
            }

            stride=nstride;
            assert(n<=stride);

            /*#ifndef NDEBUG
            for(unsigned i=0; i<n; i++){
                auto g=get_bead(i);
                if(memcmp(&g, &prev[i], sizeof(packed_bead))){
                    assert(0);
                }
            }
            #endif*/

            #endif
           }

           packed_bead erase(unsigned index)
           {
            float *fdata=(float*)data;

            assert(index < n);
            packed_bead res;
            res.hash=data[index];
            for(int d=0; d<3; d++){
                res.x[d]=fdata[ (1+d)*stride + index ];
                res.v[d]=fdata[ (4+d)*stride + index ];
                res.f[d]=fdata[ (7+d)*stride + index ];
                assert(res.f[d]==0); // temporary
            }
            n -= 1;
            if(index != n){
                data[index] = data[n];
                for(int d=0; d<3; d++){
                    fdata[ (1+d)*stride + index ] = fdata[ (1+d)*stride + n ];
                    fdata[ (4+d)*stride + index ] = fdata[ (4+d)*stride + n ];
                    fdata[ (7+d)*stride + index ] = fdata[ (7+d)*stride + n ];
                }                
            }
            return res;
           }

           void insert(const packed_bead &b)
           {
            assert(n < MAX_BEADS);
            assert(n <= stride);
            if(n == stride){
                set_stride(stride+1);
            }
            assert(n<stride);
            float *fdata=(float*)data;
            unsigned index=n;
            assert(index < MAX_BEADS);

            data[index] = b.hash;
            for(int d=0; d<3; d++){
                fdata[ (1+d)*stride+index ] = b.x[d];
                fdata[ (4+d)*stride+index ] = b.v[d];
                fdata[ (7+d)*stride+index ] = b.f[d];
            }
            ++n;

            assert(n<=stride);
           }

           const uint32_t * get_hash_vec() const
           { return data+0; } 

           const float * get_x_vec(int d) const
           { return (const float*)data+(1+d)*stride; }

           float * get_x_vec(int d)
           { return (float*)data+(1+d)*stride; }

           const float * get_v_vec(int d) const
           { return (const float*)data+(4+d)*stride; }

           float *  get_v_vec(int d)
           { return (float*)data+(4+d)*stride; } 

           const float *  get_f_vec(int d) const
           { return (const float*)data+(7+d)*stride; } 

           float *  get_f_vec(int d)
           { return (float*)data+(7+d)*stride; } 

           packed_bead get_bead(unsigned i)
           {
            assert(i < n);
            packed_bead res;
            res.hash=get_hash_vec()[i];
            for(int d=0; d<3; d++){
                res.x[d]=get_x_vec(d)[i];
                res.v[d]=get_v_vec(d)[i];
                res.f[d]=get_f_vec(d)[i];
            }
            return res;
           }

           void newton(float origin[3], float dt, std::vector<uint8_t> &outgoing)
           {
            outgoing.clear();

            for(unsigned i=0; i<n; i++){
                auto bb=get_bead(i);
                assert(-100 < bb.f[0] && bb.f[0] < +100);
            }

            auto dtv=hn::Set(vecf_tag(), dt);
            auto dt_2v=hn::Set(vecf_tag(), dt*0.5f);

            unsigned offset=0;
            while(offset < n){
                hn::Vec<vecf_tag> xd, vd, fd;
                hn::Mask<vecf_tag> moved;

                auto valid = hn::Lt(hn::Iota(vecu_tag(), offset), hn::Set(vecu_tag(), n));
                auto invalidf = hn::RebindMask(vecf_tag(), hn::Not(valid));
                for(int d=0; d<3; d++){
                    // Debug only
                    auto hash=hn::LoadU(vecu_tag(), get_hash_vec()+offset);
                    
                    xd=hn::LoadU(vecf_tag(), get_x_vec(d)+offset);
                    vd=hn::LoadU(vecf_tag(), get_v_vec(d)+offset);
                    fd=hn::LoadU(vecf_tag(), get_f_vec(d)+offset);

                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Le( hn::Set(vecf_tag(), origin[d]),  xd ) ) ) );
                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Lt( xd, hn::Set(vecf_tag(), origin[d]+1.0f) ) ) ) );
                
                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Le( hn::Set(vecf_tag(), -100),  vd ) ) ) );
                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Lt( vd, hn::Set(vecf_tag(), +100) ) ) ) );
                    
                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Le( hn::Set(vecf_tag(), -100),  fd ) ) ) );
                    assert( hn::AllTrue(vecf_tag(), hn::Or(invalidf, hn::Lt( fd, hn::Set(vecf_tag(), +100) ) ) ) );

                    if(d==0){
                        //fprintf(stderr, "hash=%u, x[0]=%f, v[0]=%f, f[0]=%f\n", hn::GetLane(hash), hn::GetLane(xd), hn::GetLane(vd), hn::GetLane(fd));
                    }

                    auto xdo=xd, vdo=vd;

                    // update mom
                    auto half_dv = dt_2v * fd;
                    vd += half_dv;

                    // Update x
                    // auto x = b.x + b.v*dt + b.f*half(dt*dt);
                    xd += dtv * ( vd + fd*dt_2v);
                    vd += half_dv;

                    //fprintf(stderr, "         x[%d]=%f, v[%d]=%f, f[%d]=%f\n", d, hn::GetLane(xd), d, hn::GetLane(vd), d, hn::GetLane(fd));
                    

                    xd = hn::IfThenElse(invalidf, xdo, xd);
                    vd = hn::IfThenElse(invalidf, vdo, vd);

                    hn::StoreU(xd, vecf_tag(), get_x_vec(d)+offset);
                    hn::StoreU(vd, vecf_tag(), get_v_vec(d)+offset);

                    auto dim_moved = hn::Or(xd < hn::Set(vecf_tag(), origin[d]) , xd >= hn::Set(vecf_tag(), origin[d]+1.0f) );
                    if(d==0){
                        moved=dim_moved;
                    }else{
                        moved = hn::Or(moved , dim_moved);
                    }
                }
                moved=hn::And( moved, hn::RebindMask(vecf_tag(), valid) );



                assert(hn::GetLane(vd)==get_v_vec(2)[offset]);

                if(!hn::AllFalse(vecf_tag(), moved)){
                    int fnz=hn::FindFirstTrue(vecf_tag(), moved);
                    while(fnz!=-1){
                        //fprintf(stderr, "Bead left\n");
                        outgoing.push_back(offset+fnz);
                        
                        moved = hn::And( moved, hn::Not(hn::FirstN(vecf_tag(), fnz+1)) );
                        fnz=hn::FindFirstTrue(vecf_tag(), moved);
                    }
                }

                offset += K;
            }

            memset( get_f_vec(0), 0, stride*3*4);

            for(unsigned i=0; i<n; i++){
                assert(get_f_vec(0)[i]==0);
            }
           }
        };

        struct soa_packed_neighbour
        {
            intptr_t dest_shr6 : sizeof(intptr_t)*8-6 ;
            intptr_t wrap_bits : 6;

            soa_packed_neighbour(soa_packed_beads *p=0)
            {
                reset(p);
            }

            void reset(soa_packed_beads *p)
            {
                wrap_bits=0;
                dest_shr6=(intptr_t)p;

                if(get_neighbour()!=p){
                    throw std::runtime_error("Pointer doesn't round-trip.");
                }
            }

            void add_wrap_bits(unsigned bits)
            {
                wrap_bits |= bits;
            }

            soa_packed_beads *get_neighbour()
            { return (soa_packed_beads*)(dest_shr6); }

            bool has_wrap() const
            { return wrap_bits !=0; }

            int get_wrap(int d) const
            {
                int bits=(wrap_bits>>(2*d))&0x3;
                return bits==2 ? -1 : bits;
            }
        };

        struct vector_packed_beads
        {
            // Extra 1 in (NHOOD_SIZE+1) is for writing slack space
            static constexpr int MAX_NHOOD_BEADS = (NHOOD_SIZE+1) * soa_packed_beads::MAX_BEADS;

            unsigned n;

            hn::Vec<vecu_tag> _align_tag_;
            uint32_t hash[MAX_NHOOD_BEADS];
            float x[3][MAX_NHOOD_BEADS];
            float v[3][MAX_NHOOD_BEADS];
            float f[3][MAX_NHOOD_BEADS];

            LaneRotator rotator;

            vector_packed_beads()
            {
                // We need to make sure that all x values start in a valid 
                // number which is well off-grid, as we assume elsewhere that they can never be nan.
                for(int d=0; d<3; d++){
                    std::fill(x[d], x[d]+MAX_NHOOD_BEADS, -10.0);
                }
            }

            packed_bead get_bead(unsigned index)
            {
                assert(index<n);
                packed_bead res;
                res.hash=hash[index];
                for(int d=0; d<3; d++){
                    res.x[d]=x[d][index];
                    res.v[d]=v[d][index];
                    res.f[d]=f[d][index];
                }
                return res;
            }

            template<bool IsEdge>
            void gather_beads(
                float box[3],
                soa_packed_neighbour *cells
            ){
                n=0;
                for(unsigned i=0; i<NHOOD_SIZE; i++){
                    const soa_packed_beads &cell = *cells[i].get_neighbour();
                    float xdelta[3]={ 0.0f, 0.0f, 0.0f };

                    if(IsEdge && cells[i].has_wrap()){
                        for(int d=0; d<3; d++){
                            int wrap=cells[i].get_wrap(d);
                            xdelta[d] =  wrap * box[d];
                        }
                    }

                    unsigned off=0;
                    do{ // This loop body still works when cell.n==0, which is quite unlikely
                        hn::StoreU( hn::LoadU( vecu_tag(), cell.get_hash_vec()+off), vecu_tag(), hash+n );
                        #pragma GCC unroll 3
                        for(int d=0; d<3; d++){
                            auto vx = hn::LoadU( vecf_tag(), cell.get_x_vec(d)+off);
                            if(IsEdge){
                                vx = hn::Set(vecf_tag(), xdelta[d]) + vx;
                            }
                            hn::StoreU( vx, vecf_tag(), x[d]+n );
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_v_vec(d)+off), vecf_tag(), v[d]+n );
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_f_vec(d)+off), vecf_tag(), f[d]+n );
                        }
                        
                        n += std::min<unsigned>(cell.n-off, K );
                        off += K;
                    }while( off < cell.n );
                }

                // This ensures that if we read off the end of the array then 
                // dr > 1. This makes the test int the main loop cheaper.
                hn::StoreU( hn::Set(vecf_tag(), -10.0), vecf_tag(), x[0]+n);
    
                // This makes sure we don't have trailing nonsense bead types, which interfer with table lookups later
                hn::StoreU( hn::Set(vecu_tag(), 0), vecu_tag(), hash+n);
            }

            void gather_beads_v2_no_edge(
                float box[3],
                soa_packed_neighbour *cells
            ){
                unsigned acc_n=0;
                hn::Vec<vecu_tag> acc_hash;
                hn::Vec<vecf_tag> acc_x[3], acc_v[3], acc_f[3];
                unsigned dst_off=0;

                n=0;
                for(unsigned i=0; i<NHOOD_SIZE; i++){
                    const soa_packed_beads &cell = *cells[i].get_neighbour();

                    assert(!cells[i].has_wrap());
                    assert(cells[i].get_neighbour()->n <= K);

                    unsigned src_off=0;
                    hn::Vec<vecu_tag> src_hash;
                    hn::Vec<vecf_tag> src_x[3], src_v[3], src_f[3];
                    unsigned src_n=cell.n;
                
                    assert(acc_n < K);
                    // The rotation is constant as acc_n == n%K, and src_n == cell.n %K
                    auto rot=rotator.get_rot(acc_n);

                continue_src:
                    src_hash=hn::LoadU( vecu_tag(), cell.get_hash_vec()+src_off);
                    for(int d=0; d<3; d++){
                        src_x[d] = hn::LoadU( vecf_tag(), cell.get_x_vec(d)+src_off);
                        src_v[d] = hn::LoadU( vecf_tag(), cell.get_v_vec(d)+src_off);
                        src_f[d] = hn::LoadU( vecf_tag(), cell.get_f_vec(d)+src_off);
                    }
                    auto keep_mask = hn::FirstN(vecf_tag(), acc_n);

                    acc_n += src_n;
                    acc_hash=hn::IfThenElse( hn::RebindMask(vecu_tag(), keep_mask), acc_hash, src_hash );
                    for(int d=0; d<3; d++){
                        acc_x[d] = hn::IfThenElse( keep_mask, acc_x[d], src_x[d]);
                        acc_v[d] = hn::IfThenElse( keep_mask, acc_v[d], src_v[d]);
                        acc_f[d] = hn::IfThenElse( keep_mask, acc_f[d], src_f[d]);
                    }
                    
                    if( (acc_n >= K) || (i==NHOOD_SIZE-1)){
                        if(acc_n < K){ // This works without the condition too
                            acc_hash = hn::IfThenElseZero( hn::FirstN(vecu_tag(), acc_n), acc_hash );
                            acc_x[0] = hn::IfThenElse( hn::FirstN(vecf_tag(), acc_n), acc_x[0], hn::Set(vecf_tag(), -10.0f));
                        }

                        hn::Store( acc_hash, vecu_tag(), hash+n);
                        acc_hash = src_hash;
                        for(int d=0; d<3; d++){
                            hn::Store( acc_x[d], vecf_tag(), x[d]+n);
                            acc_x[d] = src_x[d];
                            hn::Store( acc_v[d], vecf_tag(), v[d]+n);
                            acc_v[d] = src_v[d];
                            hn::Store( acc_f[d], vecf_tag(), f[d]+n);
                            acc_f[d] = src_f[d];
                        }
                        n += K;
                        acc_n -= K;
                    }

                    if(src_n > K){
                        src_off += K;
                        src_n -= K;
                        goto continue_src;
                    }
                }
            }



            /*
            Adds the forces in the cells.
            Assumes we have exclusive access to destinations.
            */
            void scatter_forces(
                soa_packed_neighbour *cells
            ){
                unsigned src_offset=0;
                for(unsigned i=0; i<NHOOD_SIZE; i++){
                    soa_packed_beads &cell=*cells[i].get_neighbour();
                    float *dst=cell.get_f_vec(0);

                    // Have to do in strict dimension order, as earlier dims
                    // clobber later dims.
                    for(int d=0; d<3; d++){
                        for(unsigned v=0; v<cell.n; v+=K){
                            auto f_new=hn::LoadU( vecf_tag(), f[d]+src_offset+v );

                            hn::StoreU( f_new, vecf_tag(), dst );
                            dst += std::min<unsigned>( K, cell.n-v );
                        }
                        dst += cell.stride-cell.n;
                    }
                    src_offset += cell.n;
                }
            }
        };

        static hn::Vec<vecf_tag> default_hash(
            uint64_t t_hash,
            hn::Vec<vecu_tag> a_hash,
            hn::Vec<vecu_tag> b_hash
        ){
            uint32_t a[K], b[K];
            float r[K];
            hn::StoreU(a_hash, vecu_tag(), a);
            hn::StoreU(b_hash, vecu_tag(), b);

            for(unsigned i=0; i<K; i++){
                r[i]=dpd_maths_core::default_hash(t_hash, a[i], b[i]);
                //fprintf(stderr, " v(%llu,%u,%u) -> %f\n", t_hash, a[i], b[i], r[i]);
            }

            return hn::LoadU(vecf_tag(), r);
        }

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
        template<class Config=ExecConfig>
         void filter_and_calc_force_raw(
            float rng_scale_s32_to_u, // Converts a value in [-2^31,2^31) to correctly scaled float
            uint64_t t_hash,

            const float * HWY_RESTRICT conStrength,
            const float * HWY_RESTRICT sqrtDissStrength,

            unsigned nTotal,
            unsigned nHome,  // The first nHome are active interactors

            const uint32_t * HWY_RESTRICT hash, // top 4 bits are bead type index
            const float * HWY_RESTRICT x[3],
            const float * HWY_RESTRICT v[3],
            float * HWY_RESTRICT f[3]
        ) {
            constexpr bool TSharedSqrtDissStrength=Config::SHARED_SQRT_DISS;
            constexpr int MAX_BEAD_TYPES=Config::MAX_BEAD_TYPES;

            auto get_lane=[](auto x, int lane)
            {
                decltype(hn::GetLane(x)) lanes[K];
                memcpy(lanes, &x, sizeof(x));
                return lanes[lane];
            };

            unsigned nTotalVectors = (nTotal + Lanes(vecf_tag()) - 1) / Lanes(vecf_tag());
            unsigned nTotalPadded = nTotalVectors * Lanes(vecf_tag());

            const auto sharedSqrtDissStrength=hn::Set(vecf_tag(),sqrtDissStrength[0]);

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
                    if(!TSharedSqrtDissStrength){
                        sqrtDissStrengthTable=hn::Load(vecf_tag(), sqrtDissStrengthRow);
                    }
                }

                hn::Vec<vecf_tag> f_home_acc[3]={ hn::Zero(vecf_tag()), hn::Zero(vecf_tag()), hn::Zero(vecf_tag()) };
                for(unsigned oi=0; oi<nTotalPadded; oi+=hn::Lanes(vecf_tag())){
                    hn::Vec<vecf_tag> ox[3], dx[3];
                    for(int d=0; d<3; d++){
                        ox[d]=hn::Load(vecf_tag(), x[d]+oi);
                        dx[d]=hn::Sub(hn::Set(vecf_tag(), hx[d]) , ox[d]);
                    }
                    auto dr2=hn::MulAdd(dx[0],dx[0],hn::MulAdd(dx[1],dx[1], dx[2]*dx[2]));
                    if(0 && hn::AllTrue(vecf_tag(), dr2 >= hn::Set(vecf_tag(), 1.0f))){
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
                        if(!TSharedSqrtDissStrength){
                            sqrtDissStrengthVal = hn::TableLookupLanes(sqrtDissStrengthTable, index);
                        }
                    }else{
                        conStrengthVal = hn::GatherIndex(vecf_tag(), conStrengthRow, hn::BitCast(vecs_tag(), obead_type));
                        if(!TSharedSqrtDissStrength){
                            sqrtDissStrengthVal = hn::GatherIndex(vecf_tag(), sqrtDissStrengthRow, hn::BitCast(vecs_tag(), obead_type));
                        }
                    }
                    if(TSharedSqrtDissStrength){
                        sqrtDissStrengthVal=sharedSqrtDissStrength;
                    }
                    
                    auto conForce = conStrengthVal*wr;

                    // dx is unnormalised
                    auto rdotv_mul_dr = hn::MulAdd(dx[0],dv[0],hn::MulAdd(dx[1],dv[1], dx[2]*dv[2]));

                    auto sqrt_gammap = sqrtDissStrengthVal*wr;
                    auto negDissForce = sqrt_gammap*sqrt_gammap*rdotv_mul_dr*inv_dr;

                    hn::Vec<vecf_tag> u;  // Random number in [-1,+1]
                    hn::Vec<vecf_tag> us; // Scaled random number.

                    if(!Config::USE_HASH){
                        uint32_t tt_hash=t_hash;
                        auto ubits = hn::Set(vecu_tag(), hhash + tt_hash) ^ (ohash + hn::Set(vecu_tag(), tt_hash)); // TODO: fix
                        auto ubits_s = hn::BitCast(vecs_tag(), ubits);
                        u=hn::ConvertTo(vecf_tag(), ubits_s);
                        us = u * hn::Set(vecf_tag(), rng_scale_s32_to_u ); // TODO: Check range
                    }else{
                        u=default_hash(t_hash, hn::Set(vecu_tag(), hhash), ohash);
                        us = u * hn::Set(vecf_tag(), ldexp(rng_scale_s32_to_u,31));
                    }

//                    auto randForce = sqrt_gammap * u;
//                    auto scaled_force = (conForce - negDissForce + randForce ) * inv_dr;
                    auto scaled_force = (conForce + hn::MulAdd(sqrt_gammap,us, negDissForce)) * inv_dr;

                    auto valid_r2_lt_1 = hn::And(dr < hn::Set(vecf_tag(), 1.0f) , dr > hn::Set(vecf_tag(), 0.00001f));
                    hn::Vec<vecu_tag> oi_vec=hn::Iota(vecu_tag(), oi);
                    auto valid_skip_done = hn::RebindMask(vecf_tag(), hn::Set(vecu_tag(), hi) < oi_vec );
                    auto valid = hn::And(valid_r2_lt_1 , valid_skip_done);
                    
                    scaled_force = hn::IfThenElseZero( valid, scaled_force);

                    /*for(unsigned i=0; i<K; i++){
                        dr = hn::IfThenElseZero(valid, dr);
                        if(0.0 < get_lane(dr,i) && get_lane(dr,i)<1.0f){
                            assert(get_lane(scaled_force,i)!=0);
                        }
                    }*/

                    hn::Vec<vecf_tag> f_home[3];
                    for(int d=0; d<3; d++){
                        f_home[d] = scaled_force*dx[d];
                        f_home_acc[d] += f_home[d];
                        
                        float * HWY_RESTRICT ff_o_ptr = f[d]+oi;
                        auto ff_o = hn::Load(vecf_tag(), ff_o_ptr);
                        // dx is unnormalised, but scale_force incorporates inv_dr
                        auto ff_n=ff_o - f_home[d];
                        hn::Store(ff_n, vecf_tag(), ff_o_ptr);
                        //fprintf(stderr, " o_ptr=%p\n", ff_o_ptr);
                    }

                    if(ForceLogging::logger()){
                        auto logger=ForceLogging::logger();
                        auto valid_bits=hn::VecFromMask(hn::RebindMask(vecu_tag(), valid));
                        BeadHash hhash_b{hhash};
                        for(unsigned ooi=0; ooi < K; ooi++){
                            BeadHash oohash{get_lane(ohash,ooi)};
                            unsigned ooi_vec=get_lane(oi_vec,ooi);
                        
                            //assert(hhash!=134217793);

                            //assert(oohash.hash!=134217793 || hhash_b.hash!=134217788);
                            //assert(hhash_b.hash!=134217793 || oohash.hash!=134217788);


                            //std::cerr<<"Valid["<<ooi<<"]="<<get_lane(valid_bits,ooi)<<"\n";
                            if(!get_lane(valid_bits,ooi)){
                                continue;
                            }
                            //std::cerr<<"Wibble\n";

                            
                            double ddx[3]={ get_lane(dx[0],ooi),get_lane(dx[1],ooi),get_lane(dx[2],ooi)};
                           logger->LogBeadPairProperty(hhash_b,oohash,"dx", 3,ddx);
                            for(int d=0; d<3; d++){ ddx[d] = -ddx[d]; }
                           logger->LogBeadPairProperty(oohash,hhash_b,"dx", 3,ddx);
                           
                            double ddr=get_lane(dr,ooi);
                           logger->LogBeadPairProperty(hhash_b,oohash,"dr", 1,&ddr);
                            logger->LogBeadPairProperty(oohash,hhash_b,"dr", 1,&ddr);
                        
                            double ddv[3]={get_lane(dv[0],ooi),get_lane(dv[1],ooi),get_lane(dv[2],ooi)};
                           logger->LogBeadPairProperty(hhash_b,oohash,"dv", 3,ddv);

                            double dd=get_lane(sqrtDissStrengthVal, ooi);
                           logger->LogBeadPairProperty(hhash_b, oohash,"dpd-diss-strength",dd*dd);

                            //logger->LogBeadPairProperty(hhash_b,oohash,"dpd-invrootdt", get_lane(scaled_inv_root_dt));
                            double gamma_p=pow(get_lane(sqrt_gammap,ooi),2);
                            logger->LogBeadPairProperty(hhash_b,oohash,"dpd-gammap", gamma_p);
                            logger->LogBeadPairProperty(hhash_b,oohash,"dpd-rng", get_lane(u, ooi));
                            logger->LogBeadPairProperty(hhash_b,oohash,"dpd-con", get_lane(conForce,ooi));
                            logger->LogBeadPairProperty(hhash_b,oohash,"dpd-diss", -get_lane(negDissForce,ooi));
                            //logger->LogBeadPairProperty(hhash_b,oohash,"dpd-rng-scale",1, &randScale);
                            //logger->LogBeadPairProperty(hhash_b,oohash,"dpd-rand",1, &randForce);
                            double ff[3]={get_lane(f_home[0],ooi),get_lane(f_home[1],ooi),get_lane(f_home[2],ooi)};
                            logger->LogBeadPairProperty(hhash_b,oohash,"f_next_dpd", 3,ff);
                            for(int d=0; d<3; d++){ ff[d]=-ff[d]; }
                            logger->LogBeadPairProperty(oohash,hhash_b,"f_next_dpd", 3,ff);
                        }

                    }
                }

                for(int d=0; d<3; d++){
                    float fsum=hn::GetLane(hn::SumOfLanes(vecf_tag(), f_home_acc[d]));

                    f[d][hi] += fsum;
                    assert(-100 < f[d][hi] && f[d][hi] < +100);
                }
            }
        };

        template<class Config=ExecConfig>
        void filter_and_calc_force(
            float rng_scale_s32_to_u,
            uint64_t t_hash,
            float box[3],
            bool is_edge,

            const float *conStrength,
            const float *sqrtDissStrength,

            soa_packed_neighbour *cells,

            vector_packed_beads *working
        ){
            assert(cells[0].get_neighbour()->n>0); // Should avoid calling if no home beads

            if(is_edge){
                working->gather_beads<true>(box, cells);
            }else if(0){
                working->gather_beads_v2_no_edge(box, cells);
            }else{
                working->gather_beads<false>(box, cells);
            }

            const float * __restrict__ xx[3]={working->x[0],working->x[1],working->x[2]};
            //fprintf(stderr, "  p[0]=%p\n", xx[0]);
            const float * __restrict__ vv[3]={working->v[0],working->v[1],working->v[2]};
            float * __restrict__ ff[3]={working->f[0],working->f[1],working->f[2]};
            //fprintf(stderr, "  p[0]=%p\n", ff[0]);
            
            filter_and_calc_force_raw<Config>(
                rng_scale_s32_to_u,
                t_hash,
                
                conStrength,
                sqrtDissStrength,

                working->n,
                cells[0].get_neighbour()->n,
                
                working->hash,
                xx,
                vv,
                ff
            );


            working->scatter_forces(cells);

        }
    }

};
HWY_AFTER_NAMESPACE();

#endif
