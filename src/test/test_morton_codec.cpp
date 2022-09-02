#include "dpd/core/morton_codec.hpp"

#include <iostream>

morton_codec::position_t backwards(uint64_t i)
{
    morton_codec::position_t res;
    unsigned o=0;
    while(i){
        res.x |= (i&1)<<o;
        i>>=1;
        res.y |= (i&1)<<o;
        i>>=1;
        res.z |= (i&1)<<o;
        i>>=1;
        o+=1;
    }
    return res;
}

int main()
{
    morton_codec codec;

    for(uint64_t lin_ref=0; lin_ref<(1<<30); lin_ref++){
        auto pos_ref=backwards(lin_ref);

        uint64_t lin_got=codec(pos_ref);
        if(lin_got!=lin_ref){
            std::cerr<<"Forwards fail : lin_ref="<<lin_ref<<", pos_ref="<<pos_ref.x<<","<<pos_ref.y<<","<<pos_ref.z<<" lin_got="<<lin_got<<"\n";
            exit(1);
        }
        auto pos_got=codec(lin_ref);
        if(!(pos_got==pos_ref)){
            std::cerr<<"Backwards fail : lin_ref="<<lin_ref<<", pos_ref="<<pos_ref.x<<","<<pos_ref.y<<","<<pos_ref.z<<" lin_got="<<lin_got<<", pos_got="<<pos_got.x<<","<<pos_got.y<<","<<pos_got.z<<"\n";
            exit(1);
        }
    }
}