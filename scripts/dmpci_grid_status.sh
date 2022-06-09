#!/bin/bash

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

for i in ${DMPCI_FILES} ; do
    NAME=$(basename $i)
    NAME=${NAME#dmpci.}

    STATUS_FILE=${DMPCI_GRID_NAME}/working/${NAME}/${NAME}.status.txt

    if [[ ! -f ${STATUS_FILE} || -z ${STATUS_FILE} ]] ; then
        printf "%16s, %16s, Ready\n" "${DMPCI_GRID_NAME}" "${NAME}"
    else
        STATUS=$(cat ${STATUS_FILE})
        printf "%16s, %16s, %s\n" "${DMPCI_GRID_NAME}" "${NAME}" "${STATUS}"
    fi

done