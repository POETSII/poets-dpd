#!/bin/bash

set -eou pipefail

STATE_FILE="$1"

NSTEPS=300
if [[ "$#" -gt 1 ]] ; then
    NSTEPS="$2"
fi
WORKING=""
if [[ "$#" -gt 2 ]] ; then
    WORKING="$3"
fi
    
GRAPH_SCHEMA_DIR=../graph_schema

if [[ "$WORKING" == "" ]] ; then
    mkdir -p .working
    WORKING=$(mktemp -d -p .working "xml_orch.XXXXXXXXXX" )
fi
>&2 echo "WORKING=$WORKING"

{
    make bin/create_xml_v5_graph_instance bin/extract_state_from_orch_log

    >&2 echo "Creating graph instance"
    echo "bin/create_xml_v5_graph_instance $NSTEPS 1 < ${STATE_FILE} > ${WORKING}/in.xml"
    bin/create_xml_v5_graph_instance $NSTEPS 1 < ${STATE_FILE} > ${WORKING}/in.xml
    >&2 echo "...created"

    graph_type_id=$( ${GRAPH_SCHEMA_DIR}/tools/print_graph_type_id.py ${WORKING}/in.xml )

    >&2 echo "Graph type id = ${graph_type_id}"

    if [[ ! -f "providers/${graph_type_id}.graph.so" ]] ; then
        mkdir -p providers
        ${GRAPH_SCHEMA_DIR}/tools/compile_graph_as_provider.sh ${WORKING}/in.xml \
            --release \
            --working-dir ${WORKING} \
            --output-dir ${WORKING} \
            -std=c++17
        mv -b -f ${WORKING}/${graph_type_id}.graph.so providers/${graph_type_id}.graph.so
    fi

    >&2 echo "Running in graph_schema epoch_sim, log ${WORKING}/epoch_sim.out.log"
    POETS_PROVIDER_PATH=providers ${GRAPH_SCHEMA_DIR}/bin/epoch_sim --log-level 3 --stats-delta 20  ${WORKING}/in.xml 2> ${WORKING}/epoch_sim.out.log
    grep _HANDLER_EXIT_SUCCESS_9be65737_ ${WORKING}/epoch_sim.out.log || { >&2 echo " failed."; exit 1; } 
 
    for STRATEGY in LIFO FIFO Random ; do
    >&2 echo "Running in graph_schema graph_sim, strategy {STRATEGY} log ${WORKING}/graph_sim.${STRATEGY}.out.log"
        POETS_PROVIDER_PATH=providers ${GRAPH_SCHEMA_DIR}/bin/graph_sim --strategy ${STRATEGY} --log-level 3 ${WORKING}/in.xml 2> ${WORKING}/graph_sim.${STRATEGY}.out.log
        grep _HANDLER_EXIT_SUCCESS_9be65737_ ${WORKING}/graph_sim.${STRATEGY}.out.log || { >&2 echo "${STRATEGY} failed."; exit 1; } 
    done

    if [[ ! -f "providers/${graph_type_id}.poems" ]] ; then 
        mkdir -p providers 
        PYTHONPATH=${GRAPH_SCHEMA_DIR}/tools ${GRAPH_SCHEMA_DIR}/tools/poems/compile_poems_sim.sh -o providers/${graph_type_id}.poems \
            --release --working-dir ${WORKING} \
            ${WORKING}/in.xml
    fi 

    providers/${graph_type_id}.poems --log-level 4 ${WORKING}/in.xml 2> ${WORKING}/poems.out.log
    grep _HANDLER_EXIT_SUCCESS_9be65737_ ${WORKING}/poems.out.log || { >&2 echo "POEMS failed."; exit 1; }

} | tee ${WORKING}/build.log
