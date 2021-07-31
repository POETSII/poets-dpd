CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall -O0
CPPFLAGS += -Wno-unused-variable -fmax-errors=2
LDFLAGS += -fuse-ld=gold

CPPFLAGS += -O3 -march=native
#CPPFLAGS += -DNDEBUG=1
#CPPFLAGS += -fsanitize=address -fsanitize=undefined

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

bin/% : obj/%.o
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

bin/test_hash : LDLIBS += -ltestu01