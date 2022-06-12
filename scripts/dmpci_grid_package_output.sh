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

    if ls ${DMPCI_GRID_NAME}/output/${NAME}/${NAME}-.*.rst.gz > /dev/null ; then
        for i in ${DMPCI_GRID_NAME}/output/${NAME}/${NAME}-.*.rst.gz ; do
            nn=$(echo $i | sed -r -e "s/^.*[.]0*([^.]+).rst.gz$/\1/g")
            if [[ "${nn}" != "$i" ]] ; then
                mv $i ${DMPCI_GRID_NAME}/output/${NAME}/dmpccs.${NAME}.con.${nn}.rst.gz
            fi 
        done
    fi

    if ls ${DMPCI_GRID_NAME}/output/${NAME}/*.png > /dev/null ; then
        (cd ${DMPCI_GRID_NAME}/output/${NAME} && ffmpeg -y -framerate 10 -pattern_type glob -i '*.png' -r 10 ${NAME}.mp4 )
    fi
done

if ls ${DMPCI_GRID_NAME}/output/${NAME}/*.png > /dev/null ; then
    (cd ${DMPCI_GRID_DIR}/output && ${POETS_DPD_DIR}/scripts/make_2d_mosaic.py '*/*.png')

    if ls ${DMPCI_GRID_NAME}/output/out.*.png > /dev/null ; then
        (cd ${DMPCI_GRID_NAME}/output && ffmpeg -y -framerate 10 -pattern_type glob -i '*.png' -r 10 ${NAME}.mp4 )
    fi
fi


NAME="${DMPCI_GRID_NAME}-$(date +%Y-%m-%d--%H-%M-%S)"

zip ${NAME} -q -0 -r ${DMPCI_GRID_NAME}/input ${DMPCI_GRID_NAME}/output