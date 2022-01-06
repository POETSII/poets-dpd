SHELL=/bin/bash

POETS_DPD_DIR ?= ../..

NS:=$(shell seq 4 28)
SIZES:=$(foreach n,$(NS),$(n)x$(n)x$(n))
GENERATORS=stationary_water
STEPS:=1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192
STEPS:=1 8192

define experiment_state
# $1 - Generator
# $2 - WxHxD

inputs/$1_$2.state.gz :
	mkdir -p inputs
	$(POETS_DPD_DIR)/bin/create_state/create_state_$1 $(subst x, ,$2) | gzip - > $$@

all_states : inputs/$1_$2.state.gz
all_states_$1 : inputs/$1_$2.state.gz
all_states_$2 : inputs/$1_$2.state.gz

endef

$(foreach g,$(GENERATORS),$(foreach s,$(SIZES),$(eval $(call experiment_state,$(g),$(s)))))


define experiment_xml_input
# $1 - Generator
# $2 - WxHxD
# $3 - steps

inputs/$1_$2_$3.xml.gz : inputs/$1_$2.state.gz
	gunzip -k -c $$< | (cd $(POETS_DPD_DIR) && bin/create_xml_v5_graph_instance $3 ) | gzip - > $$@

all_xml_input : inputs/$1_$2_$3.xml.gz
ALL_XML_INPUT += inputs/$1_$2_$3.xml.gz

endef

$(foreach g,$(GENERATORS),$(foreach s,$(SIZES),$(foreach step,$(STEPS),$(eval $(call experiment_xml_input,$(g),$(s),$(step))))))

ORCH_TIMEOUT=240

define experiment_run_orchestrator
# $1 - Generator
# $2 - WxHxD
# $3 - steps

outputs/$1_$2_$3.orch.log : inputs/$1_$2_$3.xml.gz
	mkdir -p outputs
	gunzip -f -k $$<
	FP=$$$$(realpath inputs/$1_$2_$3.xml); \
	$(POETS_DPD_DIR)/scripts/run_app_standard_outputs.exp $$$${FP} ${ORCH_TIMEOUT} 2>&1 | tee $$@
	rm inputs/$1_$2_$3.xml

all_run_orchestrator_output : outputs/$1_$2_$3.orch.log

endef

define experiment_run_orchestrator_random
# $1 - Generator
# $2 - WxHxD
# $3 - steps

outputs/$1_$2_$3.orch-random.log : inputs/$1_$2_$3.xml.gz
	mkdir -p outputs
	gunzip -f -k $$<
	FP=$$$$(realpath inputs/$1_$2_$3.xml); \
	$(POETS_DPD_DIR)/scripts/run_app_standard_outputs.exp $$$${FP} ${ORCH_TIMEOUT} rand 2>&1 | tee $$@
	rm inputs/$1_$2_$3.xml

all_run_orchestrator_random_output : outputs/$1_$2_$3.orch-random.log

endef

define experiment_run_polite
# $1 - input base name (e.g. stationary_water_7x7x7)
# $2 - Number of iterations

outputs/$1_$2_$3.polite.csv : inputs/$1_$2.state.gz
	mkdir -p outputs
	FP=$$$$(realpath inputs/$1_$2.state.gz); \
	TT=$$$$(mktemp -d); \
	(cd ${POETS_DPD_DIR} && bin/run_world basic_dpd_engine_v5_raw_tinsel_hw <(gunzip -k -c $$$${FP}) $$$${TT}/out.state 1 $3 ) 2> outputs/$1_$2_$3.polite.stderr | tee $$@ ; \
	rm $$$${TT}/*.state $$$${TT}/*.vtk ; \
	rmdir $$$${TT}

all_run_polite_output : outputs/$1_$2_$3.polite.csv

endef

define experiment_run_polite_random
# $1 - input base name (e.g. stationary_water_7x7x7)
# $2 - Number of iterations

outputs/$1_$2_$3.polite-random.csv : inputs/$1_$2.state.gz
	mkdir -p outputs
	FP=$$$$(realpath inputs/$1_$2.state.gz); \
	TT=$$$$(mktemp -d); \
	(cd ${POETS_DPD_DIR} && POLITE_PLACER=random bin/run_world basic_dpd_engine_v5_raw_tinsel_hw <(gunzip -k -c $$$${FP}) $$$${TT}/out.state 1 $3 ) 2> outputs/$1_$2_$3.polite-random.stderr | tee $$@ ; \
	rm $$$${TT}/*.state $$$${TT}/*.vtk ; \
	rmdir $$$${TT}

all_run_polite_random_output : outputs/$1_$2_$3.polite-random.csv

endef

$(foreach g,$(GENERATORS),$(foreach s,$(SIZES),$(foreach step,$(STEPS), \
	$(eval $(call experiment_run_orchestrator,$(g),$(s),$(step)))       \
	$(eval $(call experiment_run_orchestrator_random,$(g),$(s),$(step)))       \
	$(eval $(call experiment_run_polite,$(g),$(s),$(step)))             \
	$(eval $(call experiment_run_polite_random,$(g),$(s),$(step)))             \
)))
