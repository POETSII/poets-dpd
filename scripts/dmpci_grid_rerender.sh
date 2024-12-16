#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

DMPCI_GRID_DIR=$1
NEW_DMPCI_GRID_NAME=$2
FILTER_COLOUR=$3

if [[ ! -d "${DMPCI_GRID_DIR}" ]] ; then
    >&2 echo "Input directory ${DMPCI_GRID_DIR} does not exist."
    exit 1
fi
DMPCI_GRID_NAME=$(basename ${DMPCI_GRID_DIR})

>&2 echo "DMPCI_GRID_DIR=$1"

if [[ ! -f "${DMPCI_GRID_DIR}/dmpci_grid.name" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_DIR} does not appear to be a dmpci grid directory"
    exit 1
fi

GOT_NAME=$(cat ${DMPCI_GRID_DIR}/dmpci_grid.name)
if [[ "${GOT_NAME}" != "${DMPCI_GRID_NAME}" ]] ; then
    >&2 echo "Directory ${DMPCI_GRID_DIR} is inconsistent (dmpci_grid.name does not match directory name)"
    exit 1
fi

if [[ -d "${NEW_DMPCI_GRID_NAME}" ]] ; then
    >&2 echo "Output directory ${NEW_DMPCI_GRID_NAME} already exists."
    exit 1
fi

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

set +e
DMPCI_FILES=$(cd ${DMPCI_GRID_DIR}/input && ls dmpci.*)
set -e
if [[ "${DMPCI_FILES}" == "" ]] ; then
    >&2 echo "No dmpci files in ${DMPCI_GRID_DIR}/input"
    exit 1
fi

mkdir -p "${NEW_DMPCI_GRID_NAME}"
echo "${NEW_DMPCI_GRID_NAME}" > "${NEW_DMPCI_GRID_NAME}/dmpci_grid.name"

mkdir -p "${NEW_DMPCI_GRID_NAME}/input"

for i in ${DMPCI_FILES} ; do
    NAME=$(basename $i)
    NAME=${NAME#dmpci.}
    mkdir -p "${NEW_DMPCI_GRID_NAME}/output/${NAME}"

    >&2 echo "cp ${DMPCI_GRID_DIR}/input/$i ${NEW_DMPCI_GRID_NAME}/input"
    cp ${DMPCI_GRID_DIR}/input/$i "${NEW_DMPCI_GRID_NAME}/input"

    >&2 echo "Looking in ${DMPCI_GRID_DIR}/output/${NAME}"
    for i in ${DMPCI_GRID_DIR}/output/${NAME}/dmpccs.${NAME}.con.*.pov.gz ; do
        nn=$(echo $i | sed -r -e "s/^.*[.]0*([^.]+).pov.gz$/\1/g")    
        POV_SRC="${NEW_DMPCI_GRID_NAME}/output/${NAME}/$(basename $i .pov.gz).pov"
        >&2 echo "${POV_SRC}"
        <$i gunzip -c -k | sed "/${FILTER_COLOUR}/d" > ${POV_SRC}
        ( cd ${NEW_DMPCI_GRID_NAME}/output/${NAME} && povray +WT4 -D -V +O${NAME}.$(printf "%09d" ${nn}).png $(basename ${POV_SRC}) )
        gzip ${POV_SRC}
    done
done

