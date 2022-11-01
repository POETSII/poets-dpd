#!/bin/bash

trap 'kill -INT -$pid' INT

ORCH=../Orchestrator
SRC=$(realpath $1)

cat > tmp.bat <<HERE
    load /app = "$SRC"
    tlink /app = *
    placement /spread = *
    place /dump = *
    compose /app = *
    deploy /app = *
    initialise /app = *
    run /app = *
HERE

BAT=$(realpath tmp.bat)


(cd ${ORCH} && ./orchestrate.sh /b=${BAT} ) 
