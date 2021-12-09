INPUTS_PATH=~/osprey-dpd/dt10_experiments/scaling/gen
POETS_DPD_PATH=../..

# https://stackoverflow.com/a/18137056
CURRENT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

$(info $(CURRENT_DIR))

#INPUTS=$(notdir $(wildcard $(INPUTS_PATH)/*))
INPUTS=$(foreach n,16 24 28 32 40 48 56 64 80 96,$(foreach m,full noangle nobond,3b16e5m5r0025_$(m)_$(n)x$(n)x$(n)))

ENGINES:=NaiveDPDEngineHalfStepTBB NaiveDPDEngineHalfMergeTBB
ENGINES:=$(ENGINES) BasicDPDEngineV5RawTinselHW BasicDPDEngineV5F22RawTinselHW BasicDPDEngineV7RawCacheTinselHW BasicDPDEngineV7RawTinselHW

BATCH_SIZE=1000
BATCH_COUNT=5

# $1 = input name
# $2 = engine
define template
results/$2/$1/$2-$1.out : $(INPUTS_PATH)/$1/$1.1000.state
	mkdir -p results/$2/$1
	(cd $(POETS_DPD_PATH) && bin/run_world $2 $$< $(CURRENT_DIR)/results/$2/$1/out- $(BATCH_COUNT) $(BATCH_SIZE) ) | tee $$@

all : results/$2/$1/$2-$1.out
all_$2 : results/$2/$1/$2-$1.out
all_$1 : results/$2/$1/$2-$1.out

endef

bb :
	echo $(INPUTS)

$(foreach e,$(ENGINES),$(foreach i,$(INPUTS),$(eval $(call template,$(i),$(e)))))
