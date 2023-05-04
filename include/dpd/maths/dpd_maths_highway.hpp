#ifndef dpd_maths_highway_hpp
#define dpd_maths_highway_hpp

#include <cstdint>
#include <cmath>
#include <cassert>

#include "dpd/maths/dpd_maths_core.hpp"

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
            assert(nstride >= n);
            assert(nstride <= MAX_BEADS);
            if(n==0){
                // do nothing
            }else if(nstride>stride){
                for(int i=9; i>0; i--){
                    std::copy(data+stride*i, data+stride*i+n, data+nstride*i);
                }
            }else if(nstride<stride){
                for(int i=1; i<=9; i++){
                    std::copy_backward(data+stride*i, data+stride*i+n, data+nstride*i);
                }
            }
            stride=nstride;
            assert(n<=stride);
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

            soa_packed_neighbour(soa_packed_beads *p)
                : dest_shr6(((intptr_t)p)>>6)
                , wrap_bits(0)
            {
                if(get_neighbour()!=p){
                    throw std::runtime_error("Not aligned on 64-byte boundary.");
                }
            }

            soa_packed_beads *get_neighbour()
            { return (soa_packed_beads*)(dest_shr6<<6); }

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
            // We have 13 neighbours, plus 1 more for slack.
            static constexpr int MAX_NHOOD_BEADS = 14 * soa_packed_beads::MAX_BEADS;

            unsigned n;

            hn::Vec<vecu_tag> _align_tag_;
            uint32_t hash[MAX_NHOOD_BEADS];
            float x[3][MAX_NHOOD_BEADS];
            float v[3][MAX_NHOOD_BEADS];
            float f[3][MAX_NHOOD_BEADS];

            vector_packed_beads()
            {
                // We need to make sure that all x values start in a valid 
                // number, as we assume elsewhere that they can never be nan.
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

            void gather_beads(
                float box[3],
                unsigned ncells,
                soa_packed_neighbour *cells
            ){
                unsigned offset = 0;
                
                n=0;
                for(unsigned i=0; i<ncells; i++){
                    const soa_packed_beads &cell = *cells[i].get_neighbour();
                    float xdelta[3]={ 0.0f, 0.0f, 0.0f };
                    if(cells[i].has_wrap()){
                        for(int d=0; d<3; d++){
                            int wrap=cells[i].get_wrap(d);
                            xdelta[d] =  wrap * box[d];
                        }
                    }

                    unsigned off=0;
                    while(off < cell.n){
                        hn::StoreU( hn::LoadU( vecu_tag(), cell.get_hash_vec()+off), vecu_tag(), hash+n );
                        for(int d=0; d<3; d++){
                            hn::StoreU( hn::Set(vecf_tag(), xdelta[d]) + hn::LoadU( vecf_tag(), cell.get_x_vec(d)+off), vecf_tag(), x[d]+n );
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_v_vec(d)+off), vecf_tag(), v[d]+n );
                            hn::StoreU( hn::LoadU( vecf_tag(), cell.get_f_vec(d)+off), vecf_tag(), f[d]+n );
                        }
                        
                        n += std::min<unsigned>(cell.n-off, K );
                        off += K;
                    }
                }

                // This ensures that if we read off the end of the array then 
                // dr > 1. This makes the test int the main loop cheaper.
                hn::StoreU( hn::Set(vecf_tag(), -10.0), vecf_tag(), x[0]+n);
            }

            /*
            Adds the forces in the cells.
            Assumes we have exclusive access to destinations.
            */
            void scatter_forces(
                unsigned ncells,
                soa_packed_neighbour *cells
            ){
                unsigned src_offset=0;
                for(unsigned i=0; i<ncells; i++){
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
                hn::Vec<vecu_tag> oi_vec=hn::Iota(vecu_tag(), 0);
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

                    hn::Vec<vecf_tag> u;

                    if(!Config::USE_HASH){
                        uint32_t tt_hash=t_hash;
                        auto ubits = hn::Set(vecu_tag(), hhash + tt_hash) ^ (ohash + hn::Set(vecu_tag(), tt_hash)); // TODO: fix
                        auto ubits_s = hn::BitCast(vecs_tag(), ubits);
                        u = hn::ConvertTo(vecf_tag(), ubits_s) * hn::Set(vecf_tag(), rng_scale_s32_to_u ); // TODO: Check range
                    }else{
                        u=default_hash(t_hash, hn::Set(vecu_tag(), hhash), ohash);
                        u = u * hn::Set(vecf_tag(), ldexp(rng_scale_s32_to_u,31));
                    }

//                    auto randForce = sqrt_gammap * u;
//                    auto scaled_force = (conForce - negDissForce + randForce ) * inv_dr;
                    auto scaled_force = (conForce - hn::MulAdd(sqrt_gammap,u, negDissForce)) * inv_dr;

                    auto valid_r2_lt_1 = hn::And(dr < hn::Set(vecf_tag(), 1.0f) , dr > hn::Set(vecf_tag(), 0.00001f));
                    auto valid_skip_done = hn::RebindMask(vecf_tag(), hn::Set(vecu_tag(), hi) < oi_vec );
                    auto valid = hn::And(valid_r2_lt_1 , valid_skip_done);
                    
                    scaled_force = hn::IfThenElseZero( valid, scaled_force);

                    /*for(unsigned i=0; i<K; i++){
                        dr = hn::IfThenElseZero(valid, dr);
                        if(0.0 < get_lane(dr,i) && get_lane(dr,i)<1.0f){
                            assert(get_lane(scaled_force,i)!=0);
                        }
                    }*/


                    for(int d=0; d<3; d++){
                        auto f_home = scaled_force*dx[d];
                        f_home_acc[d] += f_home;
                        
                        float * HWY_RESTRICT ff_o_ptr = f[d]+oi;
                        auto ff_o = hn::Load(vecf_tag(), ff_o_ptr);
                        // dx is unnormalised, but scale_force incorporates inv_dr
                        auto ff_n=ff_o - f_home;
                        hn::Store(ff_n, vecf_tag(), ff_o_ptr);
                        //fprintf(stderr, " o_ptr=%p\n", ff_o_ptr);
                    }
                }

                for(int d=0; d<3; d++){
                    float fsum=hn::GetLane(hn::SumOfLanes(vecf_tag(), f_home_acc[d]));

                    f[d][hi] += fsum;
                    assert(-100 < f[d][hi] && f[d][hi] < +100);
                }

                oi_vec += hn::Set(vecu_tag(), K);
            }
        };

        template<class Config=ExecConfig>
        void filter_and_calc_force(
            float rng_scale_s32_to_u,
            uint64_t t_hash,
            float box[3],

            const float *conStrength,
            const float *sqrtDissStrength,

            unsigned ncells,
            soa_packed_neighbour *cells,

            vector_packed_beads *working
        ){
            assert(cells[0].get_neighbour()->n>0); // Should avoid calling if no home beads

            working->gather_beads(box, ncells, cells);

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


            working->scatter_forces(ncells, cells);

        }
    }

};
HWY_AFTER_NAMESPACE();

#endif
