#ifndef tolerant_dpd_engine_half_step_tbb_hpp
#define tolerant_dpd_engine_half_step_tbb_hpp

#include "dpd/engines/naive/naive_dpd_engine_half_step_tbb.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"
#include "dpd/core/logging.hpp"

#include <cassert>
#include <cmath>
#include <array>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range3d.h"
#include "tbb/concurrent_vector.h"

/*
    Tolerant version has been unified with normal version
    now, we only vary the acceptable bond length.
*/
class TolerantDPDEngineHalfStepTBB
    : public NaiveDPDEngineHalfStepTBB
{
public:
   TolerantDPDEngineHalfStepTBB()
   {
       m_max_bond_length=100;
   }
};

#endif
