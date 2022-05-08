

SCALING_BOARD_CONFIGS=1_1x1 2_2x1 4_2x2 6_2x3 9_3x3 12_3x4 16_4x4 20_4x5 24_4x6 30_5x6 36_6x6 42_6x7 48_6x8

define strong_scaling_tinsel_water
# $1 : Experiment name
# $2 : Variant name
# $3 : engine name
# $4 : Placer name
# $5 : board config

results-strong/$1/$2.$3-$4-$5.csv :
	mkdir -p results-strong/$1
	( \
		cd $(POETS_PDPD_DIR) ; \
		export POLITE_BOARDS_X=$$$$(echo $5 | sed -E "s/[0-9]+_([0-9]+)x[0-9]+/\1/g"); \
		export POLITE_BOARDS_Y=$$$$(echo $5 | sed -E "s/[0-9]+_[0-9]+x([0-9]+)/\1/g"); \
		export POLITE_PLACER=$4; \
		export POLITE_CHATTY=1; \
		>&2 echo "BOARDS_X=$$$${POLITE_BOARDS_X}, BOARDS_Y=$$$${POLITE_BOARDS_Y}" ; \
		/usr/bin/timeout --foreground 600 bin/step_world $3 $(BENCHMARKS_DIR)/states_bin/$1/$2/$2.state.gz null 1000 \
	) > $$@ 2> results-strong/$1/$2.$3-$4-$5.stderr.log

all_strong_$1 : results-strong/$1/$2.$3-$4-$5.csv 

endef

$(foreach engine,hemi_dpd_engine_v1_raw_tinsel_hw basic_dpd_engine_v5_raw_tinsel_hw gals_dpd_engine_v1_raw_tinsel_hw, \
	$(foreach config,$(SCALING_BOARD_CONFIGS), \
		$(eval $(call strong_scaling_tinsel_water,water,water-22x22x22,$(engine),metis,$(config))) \
		$(eval $(call strong_scaling_tinsel_water,water,water-32x32x32,$(engine),metis,$(config))) \
		$(eval $(call strong_scaling_tinsel_water,water,water-46x46x46,$(engine),metis,$(config))) \
		$(eval $(call strong_scaling_tinsel_water,water,water-22x22x22,$(engine),scotch,$(config))) \
		$(eval $(call strong_scaling_tinsel_water,water,water-32x32x32,$(engine),scotch,$(config))) \
		$(eval $(call strong_scaling_tinsel_water,water,water-46x46x46,$(engine),scotch,$(config))) \
	) \
)

delete_empty_strong_csv :
	find results-strong/ -size 0 -name '*.csv' -print -delete