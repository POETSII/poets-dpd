

#include <cstdint>

/*
This assumes every particle has a compressed 3x10-Bit
version of their location.
So if we have 256^3 cells, that would be 50M beads
at a density of 3. Within each cell we get a granularity
of 4, 
*/
bool fast_reject(uint32_t x, uint32_t y)
{
    // Global constant. Could be put in CSR?
    const uint32_t length=?;

    // This could be compiled into hardware. Needs to be
    // chosen to be conservative
    const int32_t thresh=?;

    int16_t x0=(x>>0)&0x3FF;
    int16_t y0=(x>>10)&0x3FF;
    int16_t z0=(x>>20)&0x3FF;

    int16_t x1=(y>>0)&0x3FF;
    int16_t y1=(y>>10)&0x3FF;
    int16_t z1=(y>>20)&0x3FF;

    x1 += ((x1==0 && x0==length-1) ? length : 0) - ((x1==length-1 && x0==0) ? length : 0);
    y1 += ((y1==0 && y0==length-1) ? length : 0) - ((y1==length-1 && y0==0) ? length : 0);
    z1 += ((z1==0 && z0==length-1) ? length : 0) - ((z1==length-1 && z0==0) ? length : 0);

    int16_t dx=x1-x0;
    int16_t dy=y1-y0;
    int16_t dz=z1-z0;

    int32_t dist_sqr = dx*int32_t(dx) + dy*int32_t(dy) + dz*int32_t(dz);

    return (dist_sqr <= thresh);
}

