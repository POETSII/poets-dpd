#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

if [[ $# -eq 0 ]] ; then 
    for i in *.zip ; do 
        $POETS_DPD_DIR/scripts/dmpci_grid_init.sh $i
        $POETS_DPD_DIR/scripts/dmpci_grid_enqueue.sh $(basename $i .zip)
    done
else 
    for i in $@ ; do 
        $POETS_DPD_DIR/scripts/dmpci_grid_init.sh $i
        $POETS_DPD_DIR/scripts/dmpci_grid_enqueue.sh $(basename $i .zip)
    done 
fi
