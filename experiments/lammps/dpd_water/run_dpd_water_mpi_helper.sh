#!/bin/bash

PREFIX=$1

[[ "$PREFIX" != "" ]] || { echo "No prefix"; exit 1 ; }

mkdir -p $PREFIX

LS="10 12 15 18 22 26 32 38 46 56 68 83 100 121"
T=10000

for L in $LS ; do
	HALF_N=$(( L*L*L*3/2))
	echo $HALF_N
	cat dpd_water_template.txt | \
		sed "s/__L__/$L/g" | \
		sed "s/__HALF_N__/$HALF_N/" | \
		sed "s/__T__/$T/g" > $PREFIX/dpd_water_template_L${L}.in

	mpiexec -np $SLURM_NTASKS lmp -in $PREFIX/dpd_water_template_L${L}.in | tee $PREFIX/dpd_water_template_L${L}.out
done


