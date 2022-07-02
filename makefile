ifeq ($(DISABLE_RISCV),1)
ENABLE_RISCV = 0
endif

ENABLE_RISCV ?= 0

CPPFLAGS += -Iinclude -std=c++17 -g3 -W -Wall -O0

CPPFLAGS += -Wno-unused-variable -fmax-errors=2 -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable
CPPFLAGS += -Wno-class-memaccess -Wno-invalid-offsetof
CPPFLAGS += -fopenmp
LDFLAGS += -pthread
LDFLAGS += -fuse-ld=gold

# Hacks for Soton HPC/AMD systems
# Detection is hacky: I just assume the module shell function only exists on iridis
ifneq ($(shell module 2>&1),)
CPPFLAGS += -I/home/dbt1c21/packages/oneTBB-2019/include
#LDFLAGS += -L/home/dt10/.linuxbrew/lib
LDFLAGS += -L/home/dbt1c21/packages/oneTBB-2019/build/linux_intel64_gcc_cc11.1.0_libc2.17_kernel3.10.0_release/
CPPFLAGS += -I/usr/include
#CPPFLAGS += -I/home/dbt1c21/.linuxbrew/include/
CPPFLAGS += -I/scratch/dbt1c21/local/include/

endif

CPPFLAGS += -g

CPPFLAGS += -DNDEBUG=1 
CPPFLAGS += -O3

#CPPFLAGS += -fno-omit-frame-pointer

## Iridis cpuinfo:
## AMD: flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush          mmx fxsr sse sse2    ht        syscall nx mmxext fxsr_opt pdpe1gb rdtscp lm constant_tsc art                       rep_good nopl           nonstop_tsc extd_apicid aperfmperf eagerfpu pni pclmulqdq        monitor                        ssse3      fma cx16                    sse4_1 sse4_2 x2apic movbe popcnt                    aes xsave avx f16c rdrand lahf_lm cmp_legacy svm extapic cr8_legacy abm sse4a misalignsse 3dnowprefetch     osvw ibs skinit wdt tce topoext perfctr_core perfctr_nb bpext perfctr_l2 cpb cat_l3 cdp_l3 hw_pstate sme retpoline_amd                                    ssbd     ibrs ibpb stibp vmmcall                                       fsgsbase            bmi1     avx2 smep bmi2                  cqm     rdt_a                  rdseed adx smap clflushopt clwb sha_ni                            xsaveopt xsavec xgetbv1 cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local clzero irperf xsaveerptr arat npt lbrv svm_lock nrip_save tsc_scale vmcb_clean flushbyasid decodeassists pausefilter pfthreshold avic v_vmsave_vmload vgif umip overflow_recov succor smca
## Intel:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx                 pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc             aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm                                   abm                   3dnowprefetch epb                                                                              cat_l3 cdp_l3                             invpcid_single intel_ppin intel_pt ssbd mba ibrs ibpb stibp         tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm mpx rdt_a avx512f avx512dq rdseed adx smap clflushopt clwb        avx512cd avx512bw avx512vl xsaveopt xsavec xgetbv1 cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts hwp_epp pku ospke md_clear spec_ctrl intel_stibp flush_l1d

CPPFLAGS += -mavx2 -mfma

# Very rarely see improvement from PGO
#CPPFLAGS += -fprofile-generate
#CPPFLAGS += -fprofile-use

#-mavx512f -mprefer-vector-width=512
#CPPFLAGS += -fsanitize=address -fsanitize=undefined 
#CPPFLAGS += -fsanitize=undefined -fsanitize=thread
#CPPFLAGS += -fopt-info-vec -fopt-info-vec-missed

CPPFLAGS += -DTBB_PREVIEW_GLOBAL_CONTROL=1

TINSEL_ROOT = tinsel

GRAPH_SCHEMA_DIR=../graph_schema
CPPFLAGS += -I $(GRAPH_SCHEMA_DIR)/include

CPPFLAGS += -I $(TINSEL_ROOT)/include
CPPFLAGS += -I $(TINSEL_ROOT)/hostlink
CPPFLAGS += -I $(TINSEL_ROOT)/apps/POLite/util/POLiteSWSim/include/POLite

#ifneq ($(TBBROOT),)
#CPPFLAGS += -I $(TBBROOT)/include
#LDFLAGS += -L $(TBBROOT)/lib
#endif

CPPFLAGS += -I ~/local/include
LDFLAGS += -L ~/local/lib

ENGINE_LDLIBS += -ltbb


all : all_test_bin all_create_state_bin all_engines all_tools




TEST_BIN := bin/test/test_engine bin/test/test_engine_diff

all_test_bin : $(TEST_BIN)

CREATE_STATE_BIN := $(patsubst src/create_state/%.cpp,bin/create_state/%,$(wildcard src/create_state/*.cpp))

all_create_state_bin : $(CREATE_STATE_BIN)

$(CREATE_STATE_BIN) : LDLIBS += -ltbb


ENGINES := $(filter-out %.riscv,$(patsubst src/engines/%.cpp,%,$(wildcard src/engines/*.cpp)))

ifeq ($(ENABLE_RISCV),1)
ENGINES_RISCV := $(filter %.riscv,$(patsubst src/engines/%.cpp,%, $(filter-out src/engines/memcpy.riscv.cpp, $(wildcard src/engines/*.cpp))))
ENGINE_LDFLAGS += -L $(TINSEL_ROOT)/hostlink
ENGINE_LDLIBS += -l:hostlink.a
ENGINE_LDLIBS += -lmetis -lscotch
else
ENGINES_RISCV := 
ENGINES := $(filter-out %tinsel_hw,$(ENGINES))
endif

ENGINES := $(filter-out basic_dpd_engine_v6% ,$(ENGINES))

ALL_ENGINE_OBJS := $(foreach e,$(ENGINES),obj/engines/$(e).o)
ALL_ENGINE_RISCV := $(foreach e,$(ENGINES_RISCV),bin/engines/$(e).code.v bin/engines/$(e).data.v )
ALL_ENGINE_BATS := $(foreach e,$(ENGINES),obj/engines/$(e).bats)

all_engines : $(ALL_ENGINES_OBJS) $(ALL_ENGINE_RISCV) $(ALL_ENGINE_BATS)

test_results/%.txt : bin/%
	mkdir -p test_results/$(dir $*)
	-rm test_results/$*.failed
	$< | tee $@ || touch test_results/$*.failed

test : $(patsubst bin/%,test_results/%.txt,$(TEST_BIN))


TEST_BATS := $(wildcard src/*.bats src/*/*.bats )

test_results/%.bats.txt : src/%.bats
	mkdir -p test_results/$(dir $*)
	-rm test_results/$*.failed
	(bats -t $< | tee $@) || touch test_results/$*.failed

test_bats : $(patsubst src/%,test_results/%.txt,$(TEST_BATS))


test : test_bats 

obj/%.d: src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MT '$(patsubst src/%.cpp,obj/%.o,$<)' $< -MF $@

obj/%.o : src/%.cpp obj/%.d
	mkdir -p obj/$(dir $*)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

src/%.S : src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) -S $(CPPFLAGS) $(CXXFLAGS) $< -o $@


obj/engines/%.bats : src/engines/%.o src/engines/dpd_engine.bats.template
	sed -e "s/__ENGINE__/$*/g" src/engines/dpd_engine.bats.template > $@


RISCV_TOOLS = /local/orchestrator-common/orchestrator_dependencies_7/riscv32-compile-driver/bin

src/%.riscv.o : src/%.cpp
	$(RISCV_CXX) -c -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@

src/%.riscv : src/%.cpp
	$(RISCV_CXX) -x c++ -I include -DTINSEL -Wdouble-promotion -DNDEBUG=1 -Os  -ffast-math -march=rv32imf -static -nostdlib -fwhole-program -g $< -o $@


#############################################################
## RISCV build stuff

RV_TOOLS_DIR =/local/orchestrator-common/orchestrator_dependencies_7/riscv32-compile-driver/bin

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
RV_CFLAGS := $(RV_CFLAGS) -O3 -I  $(TINSEL_ROOT)/include
RV_LDFLAGS = -melf32lriscv -G 0 

RV_CFLAGS := $(RV_CFLAGS) -I include
RV_CFLAGS := $(RV_CFLAGS) -std=c++17 -DNDEBUG=1 -fwhole-program -Wdouble-promotion -g -ffast-math -Wno-unused-variable \
	 -Wno-unused-parameter

RV_CFLAGS += -I $(GRAPH_SCHEMA_DIR)/include

obj/engines/link.riscv.ld :
	TINSEL_ROOT=$(TINSEL_ROOT) \
    $(TINSEL_ROOT)/apps/POLite/util/genld.sh > $@ 

obj/engines/entry.riscv.o :
	$(RV_CPPC) $(RV_CFLAGS) -W -Wall -c -o $@ $(TINSEL_ROOT)/apps/POLite/util/entry.S

## Tweak flags to avoid it getting optimised out
obj/engines/memcpy.riscv.o : src/engines/memcpy.riscv.cpp
	$(RV_CPPC) $(filter-out -fvisibility=hidden -fwhole-program,$(RV_CFLAGS)) -W -Wall -c -o $@ src/engines/memcpy.riscv.cpp

obj/engines/%.riscv.o : src/engines/%.riscv.cpp
	mkdir -p  $(dir $@)
	$(RV_CPPC) $(RV_CFLAGS) -W -Wall -c -DTINSEL -o $@  src/engines/$*.riscv.cpp

bin/engines/%.riscv.elf : obj/engines/%.riscv.o obj/engines/entry.riscv.o obj/engines/memcpy.riscv.o obj/engines/link.riscv.ld
	mkdir -p $(dir $@)
	$(RV_LD) $(RV_LDFLAGS) -T obj/engines/link.riscv.ld -o $@ obj/engines/entry.riscv.o obj/engines/memcpy.riscv.o obj/engines/$*.riscv.o $(TINSEL_ROOT)/lib/lib.o

bin/engines/%.riscv.code.v :  bin/engines/%.riscv.elf
	$(TINSEL_ROOT)/bin/checkelf.sh $<
	$(RV_OBJCOPY) -O verilog --only-section=.text $< $@ 

bin/engines/%.riscv.data.v: bin/engines/%.riscv.elf
	$(RV_OBJCOPY) -O verilog --remove-section=.text \
                --set-section-flags .bss=alloc,load,contents $< $@ 


###############################################################

bin/% : obj/%.o
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter-out %.v,$^) -o $@ $(LDFLAGS) $(LDLIBS)

define register_engine_user
bin/$1 : $(ALL_ENGINE_OBJS) $(ALL_ENGINE_RISCV)

bin/$1 : LDFLAGS += $(ENGINE_LDFLAGS)

bin/$1 : LDLIBS += $(ENGINE_LDLIBS)

all_tools += bin/$1
endef

bin/test/test_hash : LDLIBS += -ltestu01
bin/test/test_xorshift64_avx2 : LDLIBS += -ltestu01

$(eval $(call register_engine_user,test/test_engine_diff))
$(eval $(call register_engine_user,test/test_engine))
$(eval $(call register_engine_user,benchmark_engine))
$(eval $(call register_engine_user,benchmark_engine_intervals))
$(eval $(call register_engine_user,run_world))
$(eval $(call register_engine_user,step_world))
$(eval $(call register_engine_user,test/engine_diff))
$(eval $(call register_engine_user,relax_world))

all_tools : bin/extract_state_from_orch_log
all_tools : bin/create_xml_v5_graph_instance

bin/create_xml_v5_graph_instance : LDFLAGS += -static
bin/test/test_non_conflict_grid : LDLIBS += -ltbb


# TODO : This deps don't work properly, some stuff is missing
-include $(ALL_ENGINE_OBJS:%.o=%.d)
-include $($(wildcard src/*.cpp):src/%.cpp=obj/%.d)
