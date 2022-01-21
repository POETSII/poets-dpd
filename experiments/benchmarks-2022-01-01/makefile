SHELL=/bin/bash

OSPREY_DPD_BIN=~/osprey-dpd/build/dpd-poets
POETS_DPD_DIR=../..

MAX_DIM ?= 100000

GENERATORS := 

define generator_register
# $1 - Name of generator
GENERATORS := $(sort $(GENERATORS) $1)
$1_DECLARED := 1
endef

define generic_osprey_helper
# $1 - Name of generator
# $2 - Name of variant
# $3 - dmpci generator arguments
# $4 - (optional) number of smoothing steps. Default is 100
# $5 - True (non-empty) if it contains hookean bonds
# $6 - True (non-empty) if it contains angle bonds

$(eval $(call generator_register,$1))

working/$1/$2/$2.0.state.gz : generators/$1/generate_dmpci.py
	>&2 echo "$2 - running osprey"
	mkdir -p working/$1/$2
	python3 $$< $(abspath working/$1/$2/$2.0.state.gz) $3 > working/$1/$2/dmpci.$2
	(cd working/$1/$2 && $(OSPREY_DPD_BIN) $2)

states/$1/$2/$2.state.gz : working/$1/$2/$2.0.state.gz
	>&2 echo "$2 - relaxing for $(if $4,$4,100) timesteps"
	mkdir -p states/$1/$2
	$(POETS_DPD_DIR)/bin/step_world tolerant_dpd_engine_half_step_tbb $$< $$@ $(if $4,$4,100)

states_bin/$1/$2/$2.state.gz : states/$1/$2/$2.state.gz
	>&2 echo "$2 - converting to binary and zeroing forces"
	mkdir -p states_bin/$1/$2
	$(POETS_DPD_DIR)/bin/world_state_to_binary $$< $$@ 1

$1_ALL_STATES += states/$1/$2/$2.state.gz

$1_ALL_STATES_BIN += states_bin/$1/$2/$2.state.gz

endef

include $(wildcard generators/*/*.mk)

states/%.vtk : states/%.state.gz
	$(POETS_DPD_DIR)/bin/world_state_to_vtk $< $@

states/%.state : states/%.state.gz
	gunzip -k $<

define generator_template

run_% : states/$1/%.state.gz
	mkdir -p working/$1/$$*
	$(POETS_DPD_DIR)/bin/run_world naive_dpd_engine <(gunzip -k -c $$<) working/$1/$$*/out- 100 10 10

$1_all_states : $$($1_ALL_STATES)

all_states : $1_all_states

$1_all_states_bin : $$($1_ALL_STATES_BIN)

all_states_bin : $1_all_states_bin

endef

$(foreach g,$(GENERATORS), $(eval $(call generator_template,$(g))) )