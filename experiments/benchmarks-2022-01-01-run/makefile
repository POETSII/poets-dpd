POETS_PDPD_DIR ?= ../..
BENCHMARKS_DIR ?= ~/benchmarks-dpd-2022-01-01

ifeq ($(wildcard $(BENCHMARKS_DIR)/states_bin),)
$(error Can't find directory states_bin in $(BENCHMARKS_DIR))
endif

#BENCHMARKS_DIR := $(realpath $(BENCHMARKS_DIR))
#POETS_PDPD_DIR := $(realpath $(POETS_PDPD_DIR))

EXPERIMENT_NAMES := $(shell cd $(BENCHMARKS_DIR)/states_bin && ls)

FILTER_OUT := water-100x100x100 water-83x83x83
FILTER_OUT += bilayer-177x177x32 bilayer-133x133x32
FILTER_OUT += membrane_plus_polymers-192x192x48 membrane_plus_polymers-144x144x48
FILTER_OUT += oilwater_angles-100x100x100 oilwater_angles-83x83x83
FILTER_OUT += oilwater_hookean-100x100x100 oilwater_hookean-83x83x83
FILTER_OUT += oilwater_monomer-100x100x100 oilwater_monomer-83x83x83

TIMEOUT := 360

BONDED := bilayer membrane_plus_polymers oilwater_angles oilwater_hookean


define variant_template
# $1 : Experiment name
# $2 : Variant name

results/$1/$2.basic_v5-metis.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=metis /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world basic_dpd_engine_v5_raw_tinsel_hw $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 ) > $$@ 2> results/$1/$2.basic_v5-metis.stderr.log

all_$1_basic_v5-metis : results/$1/$2.basic_v5-metis.csv
all_$1 : results/$1/$2.basic_v5-metis.csv
all_$2 : results/$1/$2.basic_v5-metis.csv
all_basic_v5-metis : results/$1/$2.basic_v5-metis.csv

results/$1/$2.basic_v5-scotch.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=scotch /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world basic_dpd_engine_v5_raw_tinsel_hw $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 ) > $$@ 2> results/$1/$2.basic_v5-scotch.stderr.log

all_$1_basic_v5-scotch : results/$1/$2.basic_v5-scotch.csv
all_$1 : results/$1/$2.basic_v5-scotch.csv
all_$2 : results/$1/$2.basic_v5-scotch.csv
all_basic_v5-scotch : results/$1/$2.basic_v5-scotch.csv

results/$1/$2.basic_v5-direct.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=direct /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world basic_dpd_engine_v5_raw_tinsel_hw $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 ) > $$@ 2> results/$1/$2.basic_v5-direct.stderr.log

all_$1_basic_v5-direct : results/$1/$2.basic_v5-direct.csv
all_$1 : results/$1/$2.basic_v5-direct.csv
all_$2 : results/$1/$2.basic_v5-direct.csv
all_basic_v5-direct : results/$1/$2.basic_v5-direct.csv

results/$1/$2.basic_v5-permutation.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=direct /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world basic_dpd_engine_v5_raw_tinsel_hw $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 ) > $$@ 2> results/$1/$2.basic_v5-permutation.stderr.log

all_$1_basic_v5-permutation : results/$1/$2.basic_v5-permutation.csv
all_$1 : results/$1/$2.basic_v5-permutation.csv
all_$2 : results/$1/$2.basic_v5-permutation.csv
all_basic_v5-permutation : results/$1/$2.basic_v5-permutation.csv

ifeq ($$(filter-out $(BONDED),$1),)

results/$1/$2.basic_v5-weighted.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=metis /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world basic_dpd_engine_v5_raw_weighted_tinsel_hw $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 ) > $$@ 2> results/$1/$2.basic_v5-weighted.stderr.log

all_$1_basic_v5-weighted : results/$1/$2.basic_v5-weighted.csv
all_$1 : results/$1/$2.basic_v5-weighted.csv
all_$2 : results/$1/$2.basic_v5-weighted.csv
all_basic_v5-weighted : results/$1/$2.basic_v5-weighted.csv

endif

endef

define experiment_template
# $1 : Experiment name

$1_VARIANTS := $(filter-out $(FILTER_OUT),$(shell cd $(BENCHMARKS_DIR)/states_bin/$1 && ls ))

all : all_$1
endef

$(foreach e,$(EXPERIMENT_NAMES),$(eval $(call experiment_template,$(e))))

$(foreach e,$(EXPERIMENT_NAMES),$(foreach v,$($(e)_VARIANTS),$(eval $(call variant_template,$(e),$(v)))))

# Useful to pick up any experiments that crashed or were interrupted and need to re-run
delete_empty_csv :
	find results/ -size 0 -name '*.csv' -print -delete
