CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall
CPPFLAGS += -Wno-unused-variable
LDFLAGS += -fuse-ld=gold

CPPFLAGS += -DNDEBUG=1 -O3 -march=native
#CPPFLAGS += -fsanitize=address

TEST_BIN := bin/test_naive_engine

all : $(TEST_BIN)

obj/%.o : src/%.cpp
	mkdir -p obj
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/% : obj/%.o
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
