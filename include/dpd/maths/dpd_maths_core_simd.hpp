#ifndef dpd_maths_core_simd_hpp
#define dpd_maths_core_simd_hpp

#include "dpd/maths/dpd_maths_core.hpp"

#include "dpd/core/hash.hpp"

#ifndef PDPD_TINSEL
#include <iostream>
#endif

#define SIMDPP_ARCH_X86_AVX2
#include "simdpp/simd.h"

namespace dpd_maths_core_simd
{

    using dpd_maths_core::default_hash;




};

#endif