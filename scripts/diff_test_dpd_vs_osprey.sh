#!/bin/bash 
set -euo pipefail
IFS=$'\n\t'

DMPCI_PATH=$1
OSPREY_PATH=$2
ENGINE=$3

[[ -f ${DMPCI_PATH} ]] || { >&2 echo "File ${DMPCI_PATH} does not exist"; exit 1; }

DMPCI_FILENAME=$(basename -- "$DMPCI_PATH")
DMPCI_EXT="${DMPCI_FILENAME##*.}"
DMPCI_BASE="${DMPCI_FILENAME%.*}"

>&2 echo "path=${DMPCI_PATH}, filename=${DMPCI_FILENAME}, ext=${DMPCI_EXT}, base=${DMPCI_BASE}"

[[ "${DMPCI_BASE}" == "dmpci" ]] || { >&2 echo "File ${DMPCI_PATH} does not have a dmpci basename. Got ${DMPCI_BASE}"; exit 1; }

grep -E "Command\s+TogglePDPDHashRandom" ${DMPCI_PATH} && { >&2 echo "File ${DMPCI_PATH} already has a TogglePDPDHashRandom command"; exit 1; }
grep -E "Command\s+ToggleForceLogging"  ${DMPCI_PATH} && { >&2 echo "File ${DMPCI_PATH} already has a ToggleForceLogging command"; exit 1; }

WORKING=$(mktemp -d --tmpdir=. .tmp-dpd-diff.XXXXXXXX)
>&2 echo "Temp dir=${WORKING} (delete on exit)"
#trap 'rm -rf -- "$WORKING"' EXIT

>&2 echo "Copying in dmpci file"
cp $DMPCI_PATH ${WORKING}

DMPCI_PATH="${WORKING}/${DMPCI_FILENAME}"
>&2 echo "Working path = ${DMPCI_PATH}"

>&2 echo "Turning on PDPDHashRandom"
echo "Command TogglePDPDHashRandom 1" >> ${DMPCI_PATH}

>&2 echo "Turning on force logging"
echo "Command ToggleForceLogging 1" >> ${DMPCI_PATH}

n=20
n2=$((n*2))

>&2 echo "Removing any changes to time step"
sed -i -r -e "/Command\s+SetTimeStepSize.+/d" ${DMPCI_PATH}

>&2 echo "Setting time steps to ${n2}, all other periods to ${n}"
sed -i -r -e "s/Time\s+[0-9]+/Time ${n2}/g"   ${DMPCI_PATH}
sed -i -r -e "s/Period\s+[0-9]+/Period ${n}/g"   ${DMPCI_PATH}

>&2 echo "Asking for state dump at time-step ${n}"
echo "Command ExportToPDPDWorldState ${n} dmpci.${n}.state" >> ${DMPCI_PATH}

>&2 echo "Running osprey"
(cd ${WORKING} && ${OSPREY_PATH} ${DMPCI_EXT}) | sort > ${WORKING}/osprey.forces

[[ ! -z ${WORKING}/osprey.forces ]] || { >&2 echo "Osprey didnt seem to produce any forces"; exit 1; }
[[ ! -z ${WORKING}/dmpci.1.state ]] || { >&2 echo "Osprey didnt seem to produce a state file"; exit 1; }

>&2 echo "Running run_world n=$((n+1))"
PDPD_LOG=${WORKING}/pdpd.forces.unsorted bin/run_world ${ENGINE} ${WORKING}/dmpci.${n}.state ${WORKING}/out 1 $((n+1))
sort ${WORKING}/pdpd.forces.unsorted > ${WORKING}/pdpd.forces
rm ${WORKING}/pdpd.forces.unsorted

>&2 echo "Diffing"

for i in $(seq ${n} ${n2}) ; do
    >&2 echo "Checking step ${i}"
    >&2 echo "  Extracting x_next for pdpd"
    grep -q -E "Prop,${i},[0-9]+,,,x_next," ${WORKING}/pdpd.forces> ${WORKING}/t1.txt
    >&2 echo "  Extracting x_next for osprey"
    grep -q -E "Prop,${i},[0-9]+,,,x_next," ${WORKING}/osprey.forces> ${WORKING}/t2.txt

    [[ ! -z  ${WORKING}/t1.txt ]] || { >&2 echo "No output for pdpd"; exit 1; }
    [[ ! -z  ${WORKING}/t2.txt ]] || { >&2 echo "No output for osprey"; exit 1; }

    >&2 echo "  diffing"
    diff ${WORKING}/t1.txt ${WORKING}/t2.txt | head -n 5
done

