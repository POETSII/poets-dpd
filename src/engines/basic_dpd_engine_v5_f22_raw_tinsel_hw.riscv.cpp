#include "dpd/engines/basic/basic_dpd_engine_v5_f22_raw_tinsel.hpp"

#include "POLiteHW.h"

using Thread = typename BasicDPDEngineV5F22RawTinsel<POLiteHW<>>::Thread;


int main()
{
  // Point thread structure at base of thread's heap
  Thread* thread = (Thread*) tinselHeapBaseSRAM();
  
  // Invoke interpreter
  thread->run();

  return 0;
}
