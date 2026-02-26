#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

DMPCI_GRID_INPUT_ZIP=$1
DMPCI_GRID_INPUT_ZIP=$(realpath ${DMPCI_GRID_INPUT_ZIP})

if [[ ! -f "${DMPCI_GRID_INPUT_ZIP}" ]] ; then
    >&2 echo "Input ${DMPCI_GRID_INPUT_ZIP} does not exist."
    exit 1
fi

if [[ "${DMPCI_GRID_INPUT_ZIP}" != *.zip ]] ; then
    >&2 echo "This expects the input to be a zip file."
    exit 1
fi

DMPCI_GRID_NAME=$(basename ${DMPCI_GRID_INPUT_ZIP} .zip)

if [[ -d ${DMPCI_GRID_NAME} ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_NAME} already exists. Cowardly exit."
    echo 1
fi

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

NEEDED_POETS_BINARIES="bin/run_world bin/relax_world bin/change_world"
for i in ${NEEDED_POETS_BINARIES} ; do
    if [[ ! -x ${POETS_DPD_DIR} ]] ; then
        >&2 echo "Missing poets binary ${i}, and this script is to cowardly to build it"
        >&2 echo ""
        >&2 echo "Go to Poets DPD directory ${POETS_DPD_DIR} and do:"
        >&2 echo "$ make ${NEEDED_POETS_BINARIES}"
    fi
done

mkdir -p ${DMPCI_GRID_NAME}/input

(cd ${DMPCI_GRID_NAME}/input && unzip ${DMPCI_GRID_INPUT_ZIP})

set +e
DMPCI_FILES=$(cd ${DMPCI_GRID_NAME}/input && ls dmpci.*)
set -e
if [[ "${DMPCI_FILES}" == "" ]] ; then
    >&2 echo "No dmpci files in root of directory"
    DMPCI_FILES=$(ls ${DMPCI_GRID_NAME}/input/*/dmpci.*)
    if [[ "${DMPCI_FILES}" == "" ]] ; then
        >&2 echo "No dmpci files in root of directory or sub-directory"
        exit 1
    fi
    >&2 echo "Found dmpci in sub-dir"
    cp ${DMPCI_FILES} ${DMPCI_GRID_NAME}/input

    DMPCI_FILES=$(cd ${DMPCI_GRID_NAME}/input && ls dmpci.*)
fi

echo 

for i in ${DMPCI_FILES} ; do
    NAME=$(basename $i)
    NAME=${NAME#dmpci.}

    mkdir -p ${DMPCI_GRID_NAME}/working/${NAME}
    mkdir -p ${DMPCI_GRID_NAME}/output/${NAME}

    DMPCI_SRC_FILE="${DMPCI_GRID_NAME}/input/${i}"

    >&2 echo "Generating script for ${NAME}, DMPCI_SRC_FILE=${DMPCI_SRC_FILE}"

    STATUS_FILE=${DMPCI_GRID_NAME}/working/${NAME}/${NAME}.status.txt
    echo "Ready" > ${STATUS_FILE}

    echo "${DMPCI_GRID_NAME}" > ${DMPCI_GRID_NAME}/dmpci_grid.name

    cat << EOF > ${DMPCI_GRID_NAME}/${NAME}.sh
#!/bin/bash
#SBATCH --partition=amd
#SBATCH --mem=16000
#SBATCH --job-name=${DMPCI_GRID_NAME}-${NAME}
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --time=24:00:00

if [[ ! -f ${STATUS_FILE} ]] ; then
    >&2 echo "No status file ${STATUS_FILE}"
    exit 1
fi

VAL=\$(cat ${STATUS_FILE})
if [[ "\${SLURM_JOB_ID}" == "" ]] ; then
    if [[ "\${VAL}" != "Ready" ]] ; then
        >&2 echo "Job ${DMPCI_GRID_NAME}-${NAME} appears to already have started."
        >&2 echo "${NAME}.status.txt:"
        >&2 cat ${STATUS_FILE}
        exit 1
    fi
else
    if [[ "\$VAL" == Queued,SlurmJobId=\${SLURM_JOB_ID},* ]] ; then
        >&2 echo "Starting job."
    else
        >&2 echo "Job ${DMPCI_GRID_NAME}-${NAME} is starting in a different job than it was queued as."
        >&2 echo "Job ID=\${SLURM_JOB_ID}".
        >&2 echo "${NAME}.status.txt:"
        >&2 cat ${STATUS_FILE}
    fi
fi

if [[ "\${SLURM_JOB_ID}" != "" ]] ; then
    echo "Running,SlurmJobId=\${SLURM_JOB_ID},\$(date -Iseconds)" > ${STATUS_FILE}
else
    echo "Running,PID=\$\$,\$(date -Iseconds)" > ${STATUS_FILE}
fi

if module list 2> /dev/null > /dev/null ; then
    # We are probably in iridis
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbt1c21/packages/oneTBB-2019/build/linux_intel64_gcc_cc11.1.0_libc2.17_kernel3.10.0_release
    module load gcc
fi

${POETS_DPD_DIR}/scripts/run_dmpci_in_poets_dpd_v2.sh $(realpath ${DMPCI_SRC_FILE} ) \
    $(realpath ${DMPCI_GRID_NAME}/working/${NAME}) \
    $(realpath ${DMPCI_GRID_NAME}/output/${NAME})
RES=\$?

if [[ "\$RES" == 0 ]] ; then
    echo "Success,\$(date -Iseconds)" > ${STATUS_FILE}
else
    echo "Failed,\$(date -Iseconds)" > ${STATUS_FILE}
fi

EOF

chmod u+x "${DMPCI_GRID_NAME}/${NAME}.sh"

done