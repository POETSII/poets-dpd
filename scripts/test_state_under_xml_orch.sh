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

    bin/create_xml_v5_graph_instance $NSTEPS < ${STATE_FILE} > ${WORKING}/in.xml

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

    POETS_PROVIDER_PATH=providers ${GRAPH_SCHEMA_DIR}/bin/epoch_sim --log-level 3 --stats-delta 20  ${WORKING}/in.xml 2> ${WORKING}/epoch_sim.out.log
    bin/extract_state_from_orch_log ${STATE_FILE} < ${WORKING}/epoch_sim.out.log > ${WORKING}/epoch_sim.out.state

    for STRATEGY in LIFO FIFO Random ; do
        POETS_PROVIDER_PATH=providers ${GRAPH_SCHEMA_DIR}/bin/graph_sim --strategy ${STRATEGY} --log-level 3 ${WORKING}/in.xml 2> ${WORKING}/graph_sim.${STRATEGY}.out.log
        bin/extract_state_from_orch_log ${STATE_FILE} < ${WORKING}/graph_sim.${STRATEGY}.out.log > ${WORKING}/grap_sim.${STRATEGY}.out.state
    done

    if [[ ! -f "providers/${graph_type_id}.poems" ]] ; then 
        mkdir -p providers 
        PYTHONPATH=${GRAPH_SCHEMA_DIR}/tools ${GRAPH_SCHEMA_DIR}/tools/poems/compile_poems_sim.sh -o providers/${graph_type_id}.poems \
            --release --working-dir ${WORKING} \
            ${WORKING}/in.xml
    fi 

    providers/${graph_type_id}.poems --log-level 4 ${WORKING}/in.xml 2> ${WORKING}/poems.out.log
    bin/extract_state_from_orch_log ${STATE_FILE} < ${WORKING}/poems.out.log > ${WORKING}/poems.out.state

} | tee ${WORKING}/build.log
