#include "dpd/engines/basic/basic_dpd_engine_v7_raw_tinsel.hpp"

#include "POLiteHW.h"

constexpr bool USE_X_CACHE=false;
using Thread = typename BasicDPDEngineV7RawTinsel<POLiteHW<>,USE_X_CACHE>::Thread;


int main()
{
  // Point thread structure at base of thread's heap
  Thread* thread = (Thread*) tinselHeapBaseSRAM();
  
  // Invoke interpreter
  thread->run();

  return 0;
}
