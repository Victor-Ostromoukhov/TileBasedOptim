#!/bin/bash

if [ $# -ne 3 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

repetition=$1
nbthreads=$2
from=$3
to=$4

#nIterations=16
#nItegrandsPerIteration=262144

nIterations=256
nItegrandsPerIteration=4096

InputDir="../../Data/Input/Tiles/"
NextIterDir="../../Data/Input/TilesLimited/"
OutputDir="../../Data/Output/Limited/SoftEllipses/OutputTiles/Repetition_${repetition}/"
TracesDir="../../Data/Output/Limited/SoftEllipses/Traces/Repetition_${repetition}/"

#~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_SoftEllipses/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 2 -g $nItegrandsPerIteration > ../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/Trace_Level_${level}.dat

 if [ $to -eq 3 ] ;
then
    infname=${InputDir}2D_0m2net_set_1_level_729.dat
    outfname=${OutputDir}2D_0m2net_set_1_level_Opt3.dat
    outputNextStep=${NextIterDir}2D_0m2net_set_1_level_OptNext3.dat
else
    infname=${NextIterDir}2D_0m2net_set_1_level_OptNext$from.dat"
    outfname={OutputDir}2D_0m2net_set_1_level_Opt$to.dat
    outputNextStep=${NextIterDir}2D_0m2net_set_1_level_OptNext$to.dat
fi


~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i infname -o $outfname --outputNextStep  --nbPoints ${lst[$level]} --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit | tee "${TracesDir}Trace_Level_${level}.dat"
