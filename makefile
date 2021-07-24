CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall -O0
CPPFLAGS += -Wno-unused-variable
LDFLAGS += -fuse-ld=gold

CPPFLAGS += -DNDEBUG=1 -O3 -march=native
#CPPFLAGS += -fsanitize=address -fsanitize=undefined

TEST_BIN := bin/test_naive_engine \
	bin/test_naive_engine_core bin/test_naive_engine_core_diff \
	bin/test_naive_engine_half_step bin/test_naive_engine_half_step_diff \
	bin/test_basic_engine \
	bin/test_basic_engine_v2


all : $(TEST_BIN)

test_results/%.txt : bin/%
	mkdir -p test_results
	-rm test_results/$*.failed
	$< | tee $@ || touch test_results/$*.failed

test : $(patsubst bin/%,test_results/%.txt,$(TEST_BIN))


obj/%.o : src/%.cpp
	mkdir -p obj
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/% : obj/%.o
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
