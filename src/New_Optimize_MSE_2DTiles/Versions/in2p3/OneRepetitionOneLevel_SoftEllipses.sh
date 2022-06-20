#!/bin/bash

if [ $# -le 4 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

repetition=$1
nbthreads=$2
from=$3
to=$4
limit=$5

#nIterations=16
#nItegrandsPerIteration=262144

nIterations=256
#nItegrandsPerIteration=4096
nItegrandsPerIteration=16384


InputDir="../../Data/Input/Tiles/"
NextIterDir="../../Data/Input/TilesLimited/SoftEllipses/OutputTiles/Repetition_${repetition}/"
OutputDir="../../Data/Output/Limited/SoftEllipses/OutputTiles/Repetition_${repetition}/"
TracesDir="../../Data/Output/Limited/SoftEllipses/Traces/Repetition_${repetition}/"
mkdir -p ${NextIterDir}
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
  
 if [ $to -eq 3 ] ;
then
    infname=${InputDir}2D_0m2net_set_1_level_729.dat
    outfname=${OutputDir}2D_0m2net_set_1_level_Opt3.dat
    outputNextStep=${NextIterDir}2D_0m2net_set_1_level_OptNext3.dat
else
    infname=${NextIterDir}2D_0m2net_set_1_level_OptNext$from.dat
    outfname=${OutputDir}2D_0m2net_set_1_level_Opt$to.dat
    outputNextStep=${NextIterDir}2D_0m2net_set_1_level_OptNext$to.dat
fi

~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outputNextStep --nbPoints $to --integrandType 2 -g $nItegrandsPerIteration --limit $limit > "${TracesDir}Trace_Level_${to}.dat"
