CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall -O0
CPPFLAGS += -Wno-unused-variable -fmax-errors=2
CPPFLAGS += -fopenmp
#LDFLAGS += -fuse-ld=gold -pthread

CPPFLAGS += -O3 -march=native
#CPPFLAGS += -DNDEBUG=1
#CPPFLAGS += -fsanitize=address -fsanitize=undefined

POLITE_DIR = /mnt/c/UserData/dt10/external/POETS/tinsel

CPPFLAGS += -I $(POLITE_DIR)/include
CPPFLAGS += -I $(POLITE_DIR)/HostLink
CPPFLAGS += -I $(POLITE_DIR)/apps/POLite/util/POLiteSWSim/include/POLite

LDFLAGS += -L $(POLITE_DIR)/hostlink
LDLIBS += -l:hostlink.a  -lmetis


TEST_BIN := bin/test_naive_engine \
	bin/test_naive_engine_core bin/test_naive_engine_core_diff \
	bin/test_naive_engine_half_step bin/test_naive_engine_half_step_diff \
	bin/test_basic_engine bin/test_basic_engine_diff  \
	bin/test_basic_engine_v2  bin/test_basic_engine_v2_diff \
	bin/test_basic_engine_v3  bin/test_basic_engine_v3_diff \


all : $(TEST_BIN)

test_results/%.txt : bin/%
	mkdir -p test_results
	-rm test_results/$*.failed
	$< | tee $@ || touch test_results/$*.failed

test : $(patsubst bin/%,test_results/%.txt,$(TEST_BIN))


obj/%.o : src/%.cpp
	mkdir -p obj
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

src/%.S : src/%.cpp
	mkdir -p obj
	$(CXX) -S $(CPPFLAGS) $(CXXFLAGS) $< -o $@

RISCV_CXX = ../orchestrator_dependencies_7/riscv32-compile-driver/bin/riscv32-unknown-elf-g++

src/%.riscv.o : src/%.cpp
	$(RISCV_CXX) -c -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@


src/%.riscv : src/%.cpp
	$(RISCV_CXX) -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@

bin/% : obj/%.o
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

bin/test_hash : LDLIBS += -ltestu01