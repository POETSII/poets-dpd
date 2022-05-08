#ifndef edge_wrap_hpp
#define edge_wrap_hpp

#include "dpd/core/vec3.hpp"

enum EdgePositionBits
{
    XLeft = 1,
    XRight = 2,
    YLeft = 4,
    YRight = 8,
    ZLeft = 16,
    ZRight = 32
};

inline uint32_t create_wrap_bits(const int32_t box[3], const int32_t loc[3])
{
    uint32_t res=0;
    for(int d=0; d<3; d++){
        if(loc[d]==0){
            res |= 1<<(2*d);
        }else if(box[d]-1==loc[d]){
            res |= 2<<(2*d);
        }
    }
    return res;
}

// Apply wrapping to a neighbouring position, based on bits for current position
inline void do_neighbour_wrap(float x[3], uint32_t bits, const float box[3])
{
    if(bits){
        const float TWO=2.0f;
        for(int d=0; d<3; d++){
            if(bits&1){
                if(x[d] >= TWO){
                    x[d] -= box[d];
                }
            }else if(bits&2){
                if(x[d] < TWO){
                    x[d] += box[d];
                }
            }
            bits >>= 2;
        }
    }
}

// Apply wrapping to a bead moving here (with bits for here), applying wrapping as needed
inline void do_movement_wrap(float x[3], uint32_t bits, const float box[3])
{
    if(bits){
        const float TWO=2.0f;
        for(int d=0; d<3; d++){
            if(bits&1){
                if(x[d] < 0){
                    x[d] += box[d];
                    if(x[d]==box[d]){
                        x[d]=0.0f; // Avoid situation where x[d]=0-eps, so x'[d]=0-eps+box[d]==box[d]
                    }
                }
            }else if(bits&2){
                if(x[d] >= box[d]){
                    x[d] -= box[d];
                }
            }
        }
    }
}

#endif
