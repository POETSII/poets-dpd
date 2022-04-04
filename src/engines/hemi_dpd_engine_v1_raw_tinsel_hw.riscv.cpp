#include "dpd/engines/hemi/hemi_dpd_engine_v1_raw_tinsel.hpp"

#include "POLiteHW.h"

using Thread = typename HemiDPDEngineV1RawTinsel<POLiteHW<>>::Thread;


int main()
{
  // Point thread structure at base of thread's heap
  Thread* thread = (Thread*) tinselHeapBaseSRAM();
  
  // Invoke interpreter
  thread->run();

  return 0;
}
