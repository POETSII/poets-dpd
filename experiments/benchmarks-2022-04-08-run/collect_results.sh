#!/bin/bash

echo "File,Model,Engine,Placement,Cube,Steps,Volume,Beads, tLoad,tCompile,tAquire,tConfigure,tExecToFirst,tExecToLast,tPerfCounters,  bpsEndToEnd,bpsCompileConfigExecExport,bpsExec"


find results/ -iname '*.csv' -size +0 -print | while read SourceFile ; do
    BaseName=$(basename $SourceFile .csv)

    Model="${BaseName%%-*}"
    BaseName="${BaseName#*-}"

    Placement="${BaseName##*-}"
    BaseName="${BaseName%-*}"

    Engine="${BaseName##*.}"
    BaseName="${BaseName%.*}"

    Cube="${BaseName}"

    
    cat ${SourceFile} | {
        IFS= read Header
       IFS=, read -r File Engine Steps Volume Beads  tLoad tCompile tAquire tConfigure tExecToFirst tExecToLast tPerfCounters   bpsEndToEnd bpsCompileConfigExecExport bpsExec 
        echo ${File},${Model},${Engine},${Placement},${Cube},${Steps},${Volume},${Beads},${tLoad},${tCompile},${tAquire},${tConfigure},${tExecToFirst},${tExecToLast},${tPerfCounters},${bpsEndToEnd},${bpsCompileConfigExecExport},${bpsExec}
    }

    #echo "$Model, $Placement, $Engine, $Cube"

    
done
