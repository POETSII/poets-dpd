#!/bin/bash

echo "File,Model,Engine,EngineVariant,Cube,Steps,Volume,Beads, tLoad,tCompile,tAquire,tConfigure,tExecToFirst,tExecToLast,tPerfCounters,  bpsEndToEnd,bpsCompileConfigExecExport,bpsExec"


find results/ -iname '*.csv' -size +0 -print | while read SourceFile ; do
    BaseName=$(basename $SourceFile .csv)
    Model="${BaseName%%-*}"
    BaseName="${BaseName#*-}"
    EngineVariant="${BaseName#*.}"
    Cube="${BaseName%%.*}"

    cat ${SourceFile} | {
        IFS= read Header
        IFS=, read -r File Engine Steps Volume Beads  tLoad tCompile tAquire tConfigure tExecToFirst tExecToLast tPerfCounters   bpsEndToEnd bpsCompileConfigExecExport bpsExec 
        echo ${File},${Model},${Engine},${EngineVariant},${Cube},${Steps},${Volume},${Beads},${tLoad},${tCompile},${tAquire},${tConfigure},${tExecToFirst},${tExecToLast},${tPerfCounters},${bpsEndToEnd},${bpsCompileConfigExecExport},${bpsExec}

    }

    
done
