#!/bin/bash

if [ $# -le 4 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

repetition=$1
nbthreads=$2
fname=$3
from=$4
to=$5
limit=$6

#nIterations=16
#nItegrandsPerIteration=262144

nIterations=16
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
    infname=${InputDir}${fname}
    outfname=${OutputDir}${fname}
    outputNextStep=${NextIterDir}${fname}
else
    infname=${NextIterDir}${fname}
    outfname=${OutputDir}${fname}
    outputNextStep=${NextIterDir}${fname}
fi

echo ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outputNextStep --nbPoints $to --integrandType 2 -g $nItegrandsPerIteration --limit $limit
#~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outputNextStep --nbPoints $to --integrandType 2 -g $nItegrandsPerIteration --limit $limit > "${TracesDir}Trace_Level_${to}.dat"
