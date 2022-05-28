#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

OSPREY_DIR=/home/dbt1c21/projects/osprey-dpd/build
POETS_DPD_DIR=/home/dbt1c21/projects/poets-dpd

ENGINE=naive_dpd_engine_half_step_tbb

SRC_DMPCI=$1
INTERVAL_COUNT=$2
INTERVAL_SIZE=$3

[[ "${SRC_DMPCI}" != "" ]] || { >&2 echo "No src dmpci." ; exit 1 ; }
[[ "${INTERVAL_COUNT}" != "" ]] || { >&2 echo "No intervale count." ; exit 1 ; }
[[ "${INTERVAL_SIZE}" != "" ]] || { >&2 echo "No intervale size." ; exit 1 ; }

BASE_NAME=${SRC_DMPCI#dmpci.}

echo ${BASE_NAME}

WORKING_DIR=".working_${BASE_NAME}"

#####################################################
## Create initial stte
if [[ ! -f ${WORKING_DIR}/${BASE_NAME}.begin.state.gz ]] ; then

mkdir -p ${WORKING_DIR}

cat ${SRC_DMPCI} | sed \
	-e "s/^Time .*$/Time 1/g" \
	-e "s/^SamplePeriod .*$/SamplePeriod 1/g" \
	-e "s/^AnalysisPeriod .*$/AnalysisPeriod 1/g" \
	-e "s/^DensityPeriod .*$/DensityPeriod 1/g" \
	-e "s/^DisplayPeriod .*$/DisplayPeriod 1/g" \
	-e "s/^RestartPeriod .*$/RestartPeriod 1/g" \
	-e '/^Grid.*$/q' > ${WORKING_DIR}/dmpci.${BASE_NAME}

echo "Command OptimisePolymerOrderingForPDPD 1" >> ${WORKING_DIR}/dmpci.${BASE_NAME}
echo "Command ExportToPDPDWorldState 1 ${BASE_NAME}.init.state.gz" >> ${WORKING_DIR}/dmpci.${BASE_NAME}

(cd ${WORKING_DIR} && ${OSPREY_DIR}/dpd-poets ${BASE_NAME} )

STSS=$(grep SetTimeStepSize ${SRC_DMPCI});

if [[ "${STSS}" != "" ]] ; then
    STSS_TIME=$(echo "${STSS}" | sed -r "s/Command\s+SetTimeStepSize\s+([0-9]+)\s+([0-9.]+)/\1/g")
    STSS_DT=$(echo "${STSS}" | sed -r "s/Command\s+SetTimeStepSize\s+([0-9]+)\s+([0-9.]+)/\2/g")
    echo "${STSS_TIME}, ${STSS_DT}"

	# The 1.1 tolerance on bond distances is a bit arbitrary, as it should vary
	# based on the number of bond pairs and the r0
    ${POETS_DPD_DIR}/bin/relax_world tolerant_dpd_engine_half_step_tbb  ${WORKING_DIR}/${BASE_NAME}.init.state.gz \
		${WORKING_DIR}/${BASE_NAME}.tmp.state.gz ${STSS_TIME} 1.1
	>&2 echo "Relaxation done"
    
    ${POETS_DPD_DIR}/bin/change_world_dt ${WORKING_DIR}/${BASE_NAME}.tmp.state.gz ${WORKING_DIR}/${BASE_NAME}.begin.state.gz ${STSS_DT}

    rm ${WORKING_DIR}/${BASE_NAME}.tmp.state.gz
else
    cp ${WORKING_DIR}/${BASE_NAME}.init.state.gz  ${WORKING_DIR}/${BASE_NAME}.begin.state.gz
fi

fi

#############################################################
## Step the initial state

>&2 echo "${WORKING_DIR}/${BASE_NAME}.begin.state.gz" 
${POETS_DPD_DIR}/bin/run_world ${ENGINE} ${WORKING_DIR}/${BASE_NAME}.begin.state.gz ${WORKING_DIR}/${BASE_NAME}- ${INTERVAL_COUNT} ${INTERVAL_SIZE} ${INTERVAL_SIZE} --gzip-snapshot --povray-render --povray-snapshot
