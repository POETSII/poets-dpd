#!/bin/bash

MODE=$1

X=$2
Y=$3
Z=$4

TARGET=1000000000

T=$(( TARGET / ( X * Y * Z * 2 ) ))
NBEADS0=$(( X*Y*Z*3 / 2 ))
NBEADS1=$(( X*Y*Z*3 / 2 ))

SUBS="-var X ${X} -var Y ${Y} -var Z ${Z} -var T ${T} -var NBEADS0 ${NBEADS0} -var NBEADS1 ${NBEADS1}"

>&2 echo "SUBS=${SUBS}"

if [[ "$MODE" == "omp" ]] ; then

    CPUS=$(nproc)

    for n in $(seq ${CPUS} -1 1 ) ; do
        OMP_NUM_THREADS=$n ~/lammps/build/lmp ${SUBS}  < dpd_two_monomers_XxYxZ_T.lmp
    done
fi