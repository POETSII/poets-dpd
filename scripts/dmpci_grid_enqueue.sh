#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if ! which sbatch ; then
    >&2 echo "sbatch command does not exist. Is SLURM running?"
    exit 1
fi

DMPCI_GRID_DIR=$1

if [[ ! -d "${DMPCI_GRID_DIR}" ]] ; then
    >&2 echo "Input directory ${DMPCI_GRID_DIR} does not exist."
    exit 1
fi

DMPCI_GRID_NAME=$(basename ${DMPCI_GRID_DIR})

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

set +e
DMPCI_FILES=$(cd ${DMPCI_GRID_NAME}/input && ls dmpci.*)
set -e
if [[ "${DMPCI_FILES}" == "" ]] ; then
    >&2 echo "No dmpci files in ${DMPCI_GRID_NAME}/input"
    exit 1
fi

for i in ${DMPCI_FILES} ; do
    NAME=$(basename $i)
    NAME=${NAME#dmpci.}

    STATUS_FILE=${DMPCI_GRID_NAME}/working/${NAME}/${NAME}.status.txt

    if [[ ! -f ${STATUS_FILE} ]] ; then
        >&2 echo "Missing status file for task ${NAME}"
        exit 1;
    fi

    STATUS=$(cat ${STATUS_FILE})
    if [[ "$STATUS" == "Ready" ]] ; then   
        >&2 echo "Queuing up ${NAME}"
        JOB_ID=$(sbatch --parsable) ${DMPCI_GRID_NAME}/${NAME}.sh --output ${DMPCI_GRID_NAME}/output/${NAME}/${NAME}.log
        # Very slight race condition here
        echo "Queued,SlurmJobId=$JOB_ID,\$(date -Iseconds)" > ${STATUS_FILE}
    fi

    STATUS=$(cat ${STATUS_FILE})
    echo "${NAME}, ${STATUS}"

done