#include "dpd/engines/basic/basic_dpd_engine_v5_raw_tinsel.hpp"

#include "POLiteHW.h"

using Engine = BasicDPDEngineV5RawNoBondsTinsel<POLiteHW<>>;
using Thread = typename Engine::Thread;

int main()
{
  // Point thread structure at base of thread's heap
  Thread* thread = (Thread*) tinselHeapBaseSRAM();
  
  // Invoke interpreter
  thread->run();

  return 0;
}
