#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

OSPREY_PATH=${OSPREY_PATH:-}
if [[ "${OSPREY_PATH}" == "" ]] ; then
	set +e
	OSPREY_PATH=$(which dpd-poets 2> /dev/null)
	RES=$?
	set -e
fi

OSPREY_SEARCH_DIRS="$HOME/projects/osprey-dpd/build $HOME/osprey-dpd/build $HOME/POETS/osprey-dpd/build"

IFS=" "
for sd in ${OSPREY_SEARCH_DIRS} ; do
	>&2 echo "Trying ${sd}/dpd-poets"
	if [[ -x "${sd}/dpd-poets" ]] ; then
		OSPREY_PATH="${sd}/dpd-poets" ;
		break
	fi
done
IFS=$'\n\t'
if [[ "${OSPREY_PATH}" == "" ]] ; then
	>&2 echo "Couldn't find POETS version of Osprey (executable dpd-poets)"
	exit 1
fi

>&2 echo "Using Osprey at ${OSPREY_PATH}"

POETS_DPD_DIR="$(dirname "$0")/.."
POETS_DPD_DIR=$(realpath ${POETS_DPD_DIR})

ENGINE=naive_dpd_engine_half_merge_tbb

SRC_DMPCI=$1
WORKING_DIR=$2
OUTPUT_DIR=$3

>&2 echo "SRC_DMPCI=${SRC_DMPCI}"
>&2 echo "WORKING_DIR=${WORKING_DIR}"
>&2 echo "OUTPUT_DIR=${OUTPUT_DIR}"

[[ "${SRC_DMPCI}" != "" ]] || { >&2 echo "No src dmpci." ; exit 1 ; }

[[ -d ${WORKING_DIR} ]] || { >&2 echo "Working directory does not exist."; exit 1 ; }
[[ -d ${OUTPUT_DIR} ]] || { >&2 echo "Output directory does not exist."; exit 1 ; }

BASE_NAME=$(basename ${SRC_DMPCI})

BASE_NAME=${BASE_NAME#dmpci.}

TOTAL_TIME=$(grep "^Time\s" ${SRC_DMPCI} | sed "s/Time//g")
SNAPSHOT_PERIOD=$(grep "^RestartPeriod" ${SRC_DMPCI} | sed "s/RestartPeriod//g")
DISPLAY_PERIOD=$(grep "^DisplayPeriod" ${SRC_DMPCI} | sed "s/DisplayPeriod//g")

>&2 echo "SNAPSHOT_PERIOD = ${SNAPSHOT_PERIOD}"

BASE_PERIOD=$(( ${SNAPSHOT_PERIOD} > ${DISPLAY_PERIOD} ? ${DISPLAY_PERIOD} : ${SNAPSHOT_PERIOD} ))

INTERVAL_COUNT=$(( ${TOTAL_TIME} / ${BASE_PERIOD} ))

>&2 echo "base name = ${BASE_NAME}"

#####################################################
## Create initial stte
if [[ ! -f ${WORKING_DIR}/${BASE_NAME}.begin.state.gz ]] ; then

	mkdir -p ${WORKING_DIR}

	cat ${SRC_DMPCI} | sed \
		-e "s/^Time\s.*$/Time 1/g" \
		-e "s/^SamplePeriod\s.*$/SamplePeriod 1/g" \
		-e "s/^AnalysisPeriod\s.*$/AnalysisPeriod 1/g" \
		-e "s/^DensityPeriod\s.*$/DensityPeriod 1/g" \
		-e "s/^DisplayPeriod\s.*$/DisplayPeriod 1/g" \
		-e "s/^RestartPeriod\s.*$/RestartPeriod 1/g" \
		-e '/^Grid.*$/q' > ${WORKING_DIR}/dmpci.${BASE_NAME}

	echo "Command OptimisePolymerOrderingForPDPD 1" >> ${WORKING_DIR}/dmpci.${BASE_NAME}
	echo "Command ExportToPDPDWorldState 1 ${BASE_NAME}.init.state.gz" >> ${WORKING_DIR}/dmpci.${BASE_NAME}

	>&2 echo "Starting osprey"
	(cd ${WORKING_DIR} && ${OSPREY_PATH} ${BASE_NAME} )
	>&2 echo "Osprey done"

	# Copy misc output files from osprey to output 
	cp ${WORKING_DIR}/*.${BASE_NAME} ${OUTPUT_DIR}

	set +e
	STSS=$(grep SetTimeStepSize ${SRC_DMPCI});
	set -e

	if [[ "${STSS}" != "" ]] ; then
		STSS_TIME=$(echo "${STSS}" | sed -r "s/Command\s+SetTimeStepSize\s+([0-9]+)\s+([0-9.]+)/\1/g")
		STSS_DT=$(echo "${STSS}" | sed -r "s/Command\s+SetTimeStepSize\s+([0-9]+)\s+([0-9.]+)/\2/g")
		>&2 echo "${STSS_TIME}, ${STSS_DT}"

		# The 1.1 tolerance on bond distances is a bit arbitrary, as it should vary
		# based on the number of bond pairs and the r0
		${POETS_DPD_DIR}/bin/relax_world tolerant_dpd_engine_half_step_tbb  ${WORKING_DIR}/${BASE_NAME}.init.state.gz \
			${WORKING_DIR}/${BASE_NAME}.tmp.state.gz ${STSS_TIME} 1.1
		>&2 echo "Relaxation done"
		
		${POETS_DPD_DIR}/bin/change_world_dt --dt ${STSS_DT} --t 0 ${WORKING_DIR}/${BASE_NAME}.tmp.state.gz ${WORKING_DIR}/${BASE_NAME}.begin.state.gz

		rm ${WORKING_DIR}/${BASE_NAME}.tmp.state.gz
	else
		${POETS_DPD_DIR}/bin/change_world_dt --t 0 ${WORKING_DIR}/${BASE_NAME}.init.state.gz  ${WORKING_DIR}/${BASE_NAME}.begin.state.gz
	fi

fi

#############################################################
## Step the initial state

>&2 echo "${WORKING_DIR}/${BASE_NAME}.begin.state.gz" 
${POETS_DPD_DIR}/bin/run_world ${ENGINE} ${WORKING_DIR}/${BASE_NAME}.begin.state.gz \
	${OUTPUT_DIR}/${BASE_NAME} \
	${INTERVAL_COUNT} ${SNAPSHOT_PERIOD} ${DISPLAY_PERIOD} \
		--gzip-snapshot --povray-render --solvent-free-snapshot

