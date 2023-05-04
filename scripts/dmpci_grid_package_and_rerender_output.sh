#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [[ $# -eq 0 ]] ; then 
    for i in $(ls */dmpci_grid.name) ; do
        DIR=$(dirname $i)
        $0 "${DIR}"
    done

    exit 0
fi

DMPCI_GRID_DIR=$1

if [[ ! -d "${DMPCI_GRID_DIR}" ]] ; then
    >&2 echo "Input directory ${DMPCI_GRID_DIR} does not exist."
    exit 1
fi

DMPCI_GRID_NAME=$(basename ${DMPCI_GRID_DIR})

if [[ ! -f "${DMPCI_GRID_DIR}/dmpci_grid.name" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_DIR} does not appear to be a dmpci grid directory"
    exit 1
fi

GOT_NAME=$(cat ${DMPCI_GRID_DIR}/dmpci_grid.name)
if [[ "${GOT_NAME}" != "${DMPCI_GRID_NAME}" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_DIR} is inconsistent (dmpci_grid.name does not match directory name)"
    exit 1
fi

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

${POETS_DPD_DIR}/scripts/dmpci_grid_package_output.sh ${DMPCI_GRID_NAME}
${POETS_DPD_DIR}/scripts/dmpci_grid_rerender.sh ${DMPCI_GRID_NAME} ${DMPCI_GRID_NAME}-no-Yellow Yellow 
${POETS_DPD_DIR}/scripts/dmpci_grid_package_output.sh ${DMPCI_GRID_NAME}-no-Yellow
