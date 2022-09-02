#MODES="omp omp_intel mpi opencl"
MODES="omp mpi kokkos"
#MODES="opencl kokkos_gpu"
#MODES="opencl_single"

VOLUMES="32 48 64 96 128 192 256"

INSTANCE_TYPE=$(ec2metadata --instance-type)
if [[ "${INSTANCE_TYPE}" == "" ]] ; then
    >&2 "Couldn't work out instance type."
    exit 1
fi

DST="out-${INSTANCE_TYPE}.csv"

rm ${DST}

for v in ${VOLUMES} ; do
    for m in ${MODES} ; do
        ./dpd_two_monomers_run_adapt.sh ${m} ${v} ${v} ${v} >> ${DST}
    done
done
