#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

DMPCI_GRID_DIR=$1

if [[ ! -d "${DMPCI_GRID_DIR}" ]] ; then
    >&2 echo "Input directory ${DMPCI_GRID_DIR} does not exist."
    exit 1
fi

DMPCI_GRID_NAME=$(basename ${DMPCI_GRID_DIR})

if [[ ! -f "${DMPCI_GRID_NAME}/dmpci_grid.name" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_NAME} does not appear to be a dmpci grid directory"
    exit 1
fi

GOT_NAME=$(cat ${DMPCI_GRID_NAME}/dmpci_grid.name)
if [[ "${GOT_NAME}" != "${DMPCI_GRID_NAME}" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_NAME} is inconsistent (dmpci_grid.name does not match directory name)"
    exit 1
fi

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

set +e
DMPCI_FILES=$(cd ${DMPCI_GRID_NAME}/input && ls dmpci.*)
set -e
if [[ "${DMPCI_FILES}" == "" ]] ; then
    >&2 echo "No dmpci files in ${DMPCI_GRID_NAME}/input"
    exit 1
fi

${POETS_DPD_DIR}/scripts/dmpci_grid_status.sh  | tee ${DMPCI_GRID_DIR}/status.txt

NAME="${DMPCI_GRID_NAME}-$(date +%Y-%m-%d--%H-%M-%S)"

zip ${NAME} -0 -r ${DMPCI_GRID_NAME}/input ${DMPCI_GRID_NAME}/output