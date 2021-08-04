CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall -O0
CPPFLAGS += -Wno-unused-variable -fmax-errors=2
CPPFLAGS += -fopenmp
#LDFLAGS += -fuse-ld=gold -pthread

CPPFLAGS += -O3 -march=native
#CPPFLAGS += -DNDEBUG=1
#CPPFLAGS += -fsanitize=address -fsanitize=undefined

TINSEL_ROOT = tinsel

CPPFLAGS += -I $(TINSEL_ROOT)/include
CPPFLAGS += -I $(TINSEL_ROOT)/hostlink
CPPFLAGS += -I $(TINSEL_ROOT)/apps/POLite/util/POLiteSWSim/include/POLite

LDFLAGS += -L $(TINSEL_ROOT)/hostlink
LDLIBS += -l:hostlink.a  -lmetis


TEST_BIN := bin/test_naive_engine \
	bin/test_naive_engine_core bin/test_naive_engine_core_diff \
	bin/test_naive_engine_half_step bin/test_naive_engine_half_step_diff \
	bin/test_basic_engine bin/test_basic_engine_diff  \
	bin/test_basic_engine_v2  bin/test_basic_engine_v2_diff \
	bin/test_basic_engine_v3  bin/test_basic_engine_v3_diff \


ENGINES := $(filter-out %.riscv,$(patsubst src/engines/%.cpp,%,$(wildcard src/engines/*.cpp)))
ENGINES_RISCV := $(filter %.riscv,$(patsubst src/engines/%.cpp,%,$(wildcard src/engines/*.cpp)))


all : $(TEST_BIN)

ALL_ENGINE_OBJS := $(foreach e,$(ENGINES),obj/engines/$(e).o)
ALL_ENGINE_RISCV := $(foreach e,$(ENGINES_RISCV),bin/engines/$(e).code.v bin/engines/$(e).data.v )

all_engines : $(ALL_ENGINES_OBJS) $(ALL_ENGINE_RISCV)

test_results/%.txt : bin/%
	mkdir -p test_results
	-rm test_results/$*.failed
	$< | tee $@ || touch test_results/$*.failed

test : $(patsubst bin/%,test_results/%.txt,$(TEST_BIN))


obj/%.o : src/%.cpp
	mkdir -p obj/$(dir $*)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

src/%.S : src/%.cpp
	mkdir -p obj
	$(CXX) -S $(CPPFLAGS) $(CXXFLAGS) $< -o $@

RISCV_TOOLS = ../orchestrator_dependencies_7/riscv32-compile-driver/bin

src/%.riscv.o : src/%.cpp
	$(RISCV_CXX) -c -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@

src/%.riscv : src/%.cpp
	$(RISCV_CXX) -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@


#############################################################
## RISCV build stuff

RV_TOOLS_DIR = /local/orchestrator-common//orchestrator_dependencies_7/riscv32-compile-driver/bin

# From tinsel/globals.mk
RV_ARCH     = rv32imf
RV_CC       = $(RV_TOOLS_DIR)/riscv32-unknown-elf-gcc
RV_CPPC     = $(RV_TOOLS_DIR)/riscv32-unknown-elf-g++
RV_LD       = $(RV_TOOLS_DIR)/riscv32-unknown-elf-ld
RV_OBJCOPY  = $(RV_TOOLS_DIR)/riscv32-unknown-elf-objcopy
RV_CFLAGS   = -mabi=ilp32 -march=$(RV_ARCH) -static -mcmodel=medany \
              -fvisibility=hidden -nostdlib -nostartfiles \
              -fsingle-precision-constant -fno-builtin-printf \
              -ffp-contract=off -fno-builtin -ffreestanding

# from POLite/utils.mk
RV_CFLAGS := $(RV_CFLAGS) -O2 -I  $(TINSEL_ROOT)/include
RV_LDFLAGS = -melf32lriscv -G 0 

RV_CFLAGS := $(RV_CFLAGS) -I include
RV_CFLAGS := $(RV_CFLAGS) -std=c++17 -DNDEBUG=1 -fwhole-program -Wdouble-promotion -g

obj/engines/link.riscv.ld :
	TINSEL_ROOT=$(TINSEL_ROOT) \
    $(TINSEL_ROOT)/apps/POLite/util/genld.sh > $@ 

obj/engines/entry.riscv.o :
	$(RV_CPPC) $(RV_CFLAGS) -W -Wall -c -o $@ $(TINSEL_ROOT)/apps/POLite/util/entry.S

obj/engines/%.riscv.o : src/engines/%.riscv.cpp
	mkdir -p obj/engines
	$(RV_CPPC) $(RV_CFLAGS) -W -Wall -c -DTINSEL -o $@  src/engines/$*.riscv.cpp

bin/engines/%.riscv.elf : obj/engines/%.riscv.o obj/engines/entry.riscv.o obj/engines/link.riscv.ld
	mkdir -p bin/engines
	$(RV_LD) $(RV_LDFLAGS) -T obj/engines/link.riscv.ld -o $@ obj/engines/entry.riscv.o obj/engines/$*.riscv.o $(TINSEL_ROOT)/lib/lib.o

bin/engines/%.riscv.code.v :  bin/engines/%.riscv.elf
	$(TINSEL_ROOT)/bin/checkelf.sh $<
	$(RV_OBJCOPY) -O verilog --only-section=.text $< $@ 

bin/engines/%.riscv.data.v: bin/engines/%.riscv.elf
	$(RV_OBJCOPY) -O verilog --remove-section=.text \
                --set-section-flags .bss=alloc,load,contents $< $@ 

###############################################################

bin/% : obj/%.o
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter-out %.v,$^) -o $@ $(LDFLAGS) $(LDLIBS)

bin/test_hash : LDLIBS += -ltestu01

bin/test_engine_diff : $(ALL_ENGINE_OBJS) $(ALL_ENGINE_RISCV)