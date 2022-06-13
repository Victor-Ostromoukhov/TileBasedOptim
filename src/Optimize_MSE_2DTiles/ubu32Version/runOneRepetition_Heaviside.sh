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

mkdir -p ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Output
mkdir -p ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Traces

nIterations=64
#nItegrandsPerIteration=16384
nItegrandsPerIteration=4096
lst=(3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
for level in ${lst[@]}
	do
        echo ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 5 -g $nItegrandsPerIteration
        ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i ../Optimize_MSE_2DTiles/Data/Optim_Data_Input/Tiles/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Output/level_${level}.dat --nbPoints ${level} --integrandType 5 -g $nItegrandsPerIteration | tee ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Traces/Trace_Level_${level}.dat
	done


#for (( level=1 ; level <= 729 ; level++ ))
