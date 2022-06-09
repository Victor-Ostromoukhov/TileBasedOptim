#!/bin/bash

if [ $# -ne 2 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads>"
	exit
fi

if [ $# -ge 1 ]; then
  repetition=$1
else
  repetition=99
fi

if [ $# -ge 2 ]; then
  nbthreads=$2
else
  nbthreads=64
fi

mkdir -p ../Repetitions_Heaviside/Repetition_${repetition}/Output
mkdir -p ../Repetitions_Heaviside/Repetition_${repetition}/Traces

nIterations=1024
#nItegrandsPerIteration=16384
nItegrandsPerIteration=4096
for (( level=47 ; level <= 729 ; level++ ))
	do
        echo ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_Heaviside/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 5 -g $nItegrandsPerIteration
        ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_Heaviside/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 5 -g $nItegrandsPerIteration | tee ../Repetitions_Heaviside/Repetition_${repetition}/Traces/Trace_Level_${level}.dat
	done
