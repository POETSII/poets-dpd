#!/bin/bash

ORCH=../Orchestrator
SRC=$(realpath $1)

cat > tmp.bat <<HERE
    load /app = "$SRC"
    tlink /app = *
    place /tfill = *
    compose /app = *

     exit /at = "end"
HERE

BAT=$(realpath tmp.bat)


(cd ${ORCH} && ./orchestrate.sh /n /b=${BAT} ) 