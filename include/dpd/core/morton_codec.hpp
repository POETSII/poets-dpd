#ifndef morton_codec_hpp
#define morton_codec_hpp

#include <array>
#include <cstdint>

/*
    Converts from 3d coordinates to morton order (Z-order)
    by interleaving bits.

    Given 4-bit co-ordinates x,y,z, we form the 12-bit morton as:

        Normal:   x3,x2,x1,x0  y3,y2,y1,y0,  z3,z2,z1,z0
        Morton:   z3,y3,z3,z2,y2,z2,z1,y1,z1,z0,y0,z0

    The converter uses two 512 entry luts to convert forwards
    and backwards. Forwards, it converts 9 bits at a time of
    each axis. Backwards it consumes 9 bits of interleaved
    order and produces 3 bits of each co-ordinate at a time.

*/
struct morton_codec
{
    struct position_t
    {
        position_t()
            : bits(0)
        {}

        position_t(uint16_t _x, uint16_t _y, uint16_t _z)
        {
            x=_x;
            y=_y;
            z=_z;
        }

        union{
            struct{
                uint16_t x, y, z, _pad;
            };
            uint64_t bits;
        };

        position_t operator|(position_t o) const
        {
            position_t res=*this;
            res.bits |= o.bits;
            return res;
        }

        position_t operator<<(unsigned dist) const
        {
            return position_t{uint16_t(x<<dist), uint16_t(y<<dist), uint16_t(z<<dist)};
        }

        bool operator==(position_t o) const
        { return bits==o.bits; }
    };

    static const unsigned LOG2K=9;
    static_assert((LOG2K%3)==0);
    
    static const unsigned K=1<<LOG2K;
    static const unsigned KMASK=K-1;

    unsigned log2n;
    std::array<uint64_t,K> lut_fwd;
    std::array<position_t,K> lut_rev;

    void make_lut_fwd()
    {
        for(uint64_t x=0; x<K; x++){
            uint64_t acc=0;
            for(unsigned i=0; i<LOG2K; i++){
                acc |= ((x>>i)&1) << (3*i);
            }
            lut_fwd[x]=acc;
        }
    }

    // Calculate for 6 bit chunk starting at bit o
    void make_lut_rev()
    {
        for(unsigned i=0; i<K; i++){
            position_t part;
            part.bits=0;
            for(unsigned j=0; j<LOG2K/3; j++){
                part.x |= ((i>>(j*3+0))&1) << j;
                part.y |= ((i>>(j*3+1))&1) << j;
                part.z |= ((i>>(j*3+2))&1) << j;
            }
            lut_rev[i]=part;
        }
    }

    morton_codec()
    {
        make_lut_fwd();
        make_lut_rev();
    }

    uint64_t operator()(position_t pos) const
    {
        return (*this)(pos.x, pos.y, pos.z);
    }

    uint64_t operator()(uint16_t x, uint16_t y, uint16_t z) const
    {
        uint64_t res=0;

        unsigned o=0;
        do{
            res |= ( (lut_fwd[x & KMASK]<<0) | (lut_fwd[y&KMASK]<<1) | (lut_fwd[z&KMASK]<<2)) << o;
            x=x>>LOG2K;
            y=y>>LOG2K;
            z=z>>LOG2K;
            o += 3*K;
        }while(x|y|z);
        return res;
    }

    position_t operator()(uint64_t x)
    {
        position_t res;
        unsigned o=0;
        while(x){
            res = res | (lut_rev[x & KMASK]<<o);
            o += LOG2K/3;
            x >>= LOG2K;
        }
        return res;
    }
};

#endif
