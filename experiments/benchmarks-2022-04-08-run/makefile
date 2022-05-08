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

TIMEOUT := 600

BONDED := bilayer membrane_plus_polymers oilwater_angles oilwater_hookean
NON_BONDED := water oilwater_monomer

define variant_engine_placer_template
# $1 : Experiment name
# $2 : Variant name
# $3 : engine name
# $4 : Placer name

results/$1/$2.$3-$4.csv :
	mkdir -p results/$1
	(cd $(POETS_PDPD_DIR) && POLITE_PLACER=$4 /usr/bin/timeout --foreground $(TIMEOUT) bin/step_world $3 $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 2000 ) > $$@ 2> results/$1/$2.$3-$4.stderr.log

all_$1 : results/$1/$2.$3-$4.csv
all_$1_$3 : results/$1/$2.$3-$4.csv
all_$1_$2 : results/$1/$2.$3-$4.csv
all_$3_$4 : results/$1/$2.$3-$4.csv

endef

define variant_template
# $1 : Experiment name
# $2 : Variant name

$(eval $(call variant_engine_placer_template,$1,$2,basic_dpd_engine_v5_raw_tinsel_hw,scotch))
$(eval $(call variant_engine_placer_template,$1,$2,basic_dpd_engine_v5_raw_tinsel_hw,metis))
$(eval $(call variant_engine_placer_template,$1,$2,basic_dpd_engine_v5_raw_tinsel_hw,direct))
$(eval $(call variant_engine_placer_template,$1,$2,basic_dpd_engine_v5_raw_tinsel_hw,permutation))

$(if $(filter-out $(NON_BONDED),$1), $(eval $(call variant_engine_placer_template,$1,$2,basic_dpd_engine_v5_raw_weighted_tinsel_hw,metis)))


$(if $(filter-out $(BONDED),$1), $(eval $(call variant_engine_placer_template,$1,$2,gals_dpd_engine_v1_raw_tinsel_hw,metis)))
$(if $(filter-out $(BONDED),$1), $(eval $(call variant_engine_placer_template,$1,$2,gals_dpd_engine_v1_raw_tinsel_hw,scotch)))
$(if $(filter-out $(BONDED),$1), $(eval $(call variant_engine_placer_template,$1,$2,hemi_dpd_engine_v1_raw_tinsel_hw,metis)))
$(if $(filter-out $(BONDED),$1), $(eval $(call variant_engine_placer_template,$1,$2,hemi_dpd_engine_v1_raw_tinsel_hw,scotch)))

#$(eval $(call variant_engine_placer_template,$1,$2,naive_dpd_engine_half_step_tbb,none))
#$(eval $(call variant_engine_placer_template,$1,$2,naive_dpd_engine_half_merge_tbb,none))

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

include strong_scaling_experiments.mk