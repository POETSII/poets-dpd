#!/bin/bash

echo "File,Model,Engine,Placement,Cube,Steps,Volume,Beads, BoardsX,BoardsY,BoardsTotal, tLoad,tCompile,tAquire,tConfigure,tExecToFirst,tExecToLast,tPerfCounters,  bpsEndToEnd,bpsCompileConfigExecExport,bpsExec"


find results-strong/ -iname '*.csv' -size +0 -print | while read SourceFile ; do
    BaseName=$(basename $SourceFile .csv)

    Model="${BaseName%%-*}"
    BaseName="${BaseName#*-}"

    Cube="${BaseName%%.*}"
    BaseName="${BaseName#*.}"

    Engine="${BaseName%%-*}"
    BaseName="${BaseName#*-}"

    Placement="${BaseName%%-*}"
    BaseName="${BaseName#*-}"

    BoardsTotal="${BaseName%%_*}"
    BaseName="${BaseName#*_}"

    BoardsX="${BaseName%%x*}"
    BaseName="${BaseName#*x}"

    BoardsY="${BaseName}"
    
    cat ${SourceFile} | {
        IFS= read Header
       IFS=, read -r File Engine Steps Volume Beads  tLoad tCompile tAquire tConfigure tExecToFirst tExecToLast tPerfCounters   bpsEndToEnd bpsCompileConfigExecExport bpsExec 
        echo ${File},${Model},${Engine},${Placement},${Cube},${Steps},${Volume},${Beads},${BoardsX},${BoardsY},${BoardsTotal},${tLoad},${tCompile},${tAquire},${tConfigure},${tExecToFirst},${tExecToLast},${tPerfCounters},${bpsEndToEnd},${bpsCompileConfigExecExport},${bpsExec}
    }

    #echo "$Model, $Placement, $Engine, $Cube"

    
done
