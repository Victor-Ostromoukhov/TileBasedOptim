#!/bin/bash

if [ $# -ne 3 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

repetition=$1
nbthreads=$2
level=$3


nIterations=1024
nItegrandsPerIteration=4096

~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_Heaviside/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 5 -g $nItegrandsPerIteration > ../Repetitions_Heaviside/Repetition_${repetition}/Traces/Trace_Level_${level}.dat
