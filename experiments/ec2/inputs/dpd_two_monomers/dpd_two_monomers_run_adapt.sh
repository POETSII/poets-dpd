#!/bin/bash

MODE=$1

X=$2
Y=$3
Z=$4

HERE=$(cd $(dirname $0) && pwd)

LMP_CPU_PATH=${HERE}/../../lammps/build_cpu/lmp
LMP_GPU_PATH=${HERE}/../../lammps/build_gpu/lmp

NBEADS0=$(( X*Y*Z*3 / 2 ))
NBEADS1=$(( X*Y*Z*3 / 2 ))

SUBS="-var X ${X} -var Y ${Y} -var Z ${Z} -var NBEADS0 ${NBEADS0} -var NBEADS1 ${NBEADS1}"

>&2 echo "SUBS=${SUBS}"

INSTANCE_TYPE=$(ec2metadata --instance-type)

mkdir logs

function find_time {
    NAME=$1
    CMD=$2

    local T=8
    local start
    local end

    while /bin/true ; do
        >&2 echo "$NAME, T=$T,  CMD=${CMD} -var T ${T} -in dpd_two_monomers_XxYxZ_T.lmp"
        start=$(date +%s)
        ${CMD} -var T ${T} -log logs/${NAME}_T${T}.log -in dpd_two_monomers_XxYxZ_T.lmp > /dev/null
        RES=$?
        end=$(date +%s)

        if [[ $RES -ne 0  ]] ; then 
            >&2 echo "Failed"
            exit 1
        fi

        ROW=$(grep "^Loop time of" logs/${NAME}_T${T}.log)
        read G_TIME G_CPUS G_STEPS G_BEADS < <(echo "$ROW" | \
            sed -r -e "s/Loop time of ([0-9.]+) on ([0-9]+) procs for ([0-9]+) steps with ([0-9]+) atoms/\1 \2 \3 \4/g" )
        
        TFLOOR=$( python3 -c "import math; print( math.floor(${G_TIME}) )" )    

        total_time=$(( end - start ))
        if [[ ${TFLOOR} -ge 10 ]] ; then 
            MBPS=$( python3 -c "print( ${G_BEADS} * ${G_STEPS} / ${G_TIME} )" )    

            echo "${INSTANCE_TYPE}, ${MODE}, ${X}, ${Y}, ${Z}, ${total_time}, ${G_TIME}, ${G_CPUS}, ${G_STEPS}, ${G_BEADS}, ${MBPS}"

            break
        fi

        >&2 echo "  -> time=$((delta))"

        T=$((T*2))
    done
}

if [[ "$MODE" == "omp" ]] ; then

    CPUS=$(nproc)
    while [[ ${CPUS} -gt 0 ]] ; do
        find_time "${MODE}_${X}_${Y}_${Z}_CPU${CPUS}" "${LMP_CPU_PATH} ${SUBS} -sf omp -pk omp ${CPUS}"
        CPUS=$((CPUS/2))
    done
elif [[ "$MODE" == "kokkos" ]] ; then
    CPUS=$(nproc)
    while [[ ${CPUS} -gt 0 ]] ; do
        find_time "${MODE}_${X}_${Y}_${Z}_CPU${CPUS}" "mpirun --use-hwthread-cpus -np ${CPUS} ${LMP_CPU_PATH} ${SUBS} -k on -sf kk"
        CPUS=$((CPUS/2))
    done
elif [[ "$MODE" == "kokkos_gpu" ]] ; then
    find_time "${MODE}_${X}_${Y}_${Z}_CPU1" "mpirun --use-hwthread-cpus -np 1 ${LMP_GPU_PATH} ${SUBS} -k on g 1 -sf kk -pk kokkos newton on neigh half  "
elif [[ "$MODE" == "intel" ]] ; then

    # NOTE: This is using half the available cores.
    # the intel implementation crashes if I try to use the full number via
    # MPI or a combination of MPI+OMP or using mode single...
    # The intel implementation crashes non-determinstically on a c6a.16xlarge,
    # and consistently if number of mpi processes is greater than 1
    export KMP_BLOCKTIME=0
    CPUS=$(( $(nproc) / 2 ))
    while [[ ${CPUS} -gt 0 ]] ; do
        find_time "${MODE}_${X}_${Y}_${Z}_CPU${CPUS}" "mpirun --use-hwthread-cpus -np ${CPUS} ${LMP_CPU_PATH} ${SUBS} -sf intel "
        CPUS=$(( CPUS / 2 ))
    done
elif [[ "$MODE" == "opt" ]] ; then

    find_time "${MODE}_${X}_${Y}_${Z}_CPU1" "${LMP_CPU_PATH} ${SUBS} -sf opt "
elif [[ "$MODE" == "mpi" ]] ; then

    CPUS=$(nproc)
    while [[ ${CPUS} -gt 0 ]] ; do
        find_time "${MODE}_${X}_${Y}_${Z}_CPU${n}" "mpirun --use-hwthread-cpus -np ${CPUS} ${LMP_CPU_PATH} ${SUBS}"
        CPUS=$((CPUS/2))
    done
elif [[ "$MODE" == "opencl" ]] ; then

    find_time "${MODE}_${X}_${Y}_${Z}_CPU${n}" "${LMP_GPU_PATH} -sf gpu ${SUBS}"
else
    >&2 echo "Didn't understand mode ${MODE}"
    exit 1
fi
