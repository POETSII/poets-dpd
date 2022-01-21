#!/bin/bash

OSPREY_DIR=$1
ENGINE=$2

[[ -d ${OSPREY_DIR} ]] || { >&2 echo "Directory ${OSPREY_DIR} does not exist. " ; exit 1; }
[[ -x ${OSPREY_DIR}/examples ]] || { >&2 echo "Directory ${OSPREY_DIR}/examples does not exist. " ; exit 1; }
[[ -x ${OSPREY_DIR}/build/dpd-poets ]] || { >&2 echo "Executable ${OSPREY_DIR}/build/dpd-poets does not exist. " ; exit 1; }

EXAMPLES="water polymer micelle"

n=0
for e in ${EXAMPLES} ; do
    >&2 echo "# Testing dmpci.$e"
    scripts/diff_test_dpd_vs_osprey.sh ${OSPREY_DIR}/examples/dmpci.${e} ${OSPREY_DIR}/build/dpd-poets ${ENGINE} 2> >(sed -E -e "s/^(.*)$/# \1/g" > /dev/stderr)
    RES=$?
    if [[ $RES -ne 0 ]] ; then 
        while IFS= read -r LINE ; do
            echo "# ${LINE}"
        done <<< "${LINE}"
        echo "not ok $n"
    else 
        echo "ok $n"
    fi

    n=$((n+1))
done
