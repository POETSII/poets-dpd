#include "fixed_math.hpp"

#include <random>
#include <iostream>

int main()
{
    ap_uint<8> x;

    x=0;
    assert(( reinterpret_as_ap_fixed<8,8>(x) == 0 ));

    x=127;
    assert(( reinterpret_as_ap_fixed<8,8>(x) == 127 ));

    x=128;
    assert(( reinterpret_as_ap_fixed<8,8>(x) == -128 ));

    x=0;
    assert(( reinterpret_as_ap_fixed<8,0>(x) == 0 ));

    x=1;
    std::cout<<reinterpret_as_ap_fixed<8,1>(x)<<"\n";
    assert(( reinterpret_as_ap_fixed<8,1>(x) == 1.0/128 ));

    x=128;
    std::cout<<reinterpret_as_ap_fixed<8,1>(x)<<"\n";
    assert(( reinterpret_as_ap_fixed<8,1>(x) == -1 ));

}