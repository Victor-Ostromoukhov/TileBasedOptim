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

mkdir -p ../Repetitions/Repetition_${repetition}/Output
mkdir -p ../Repetitions/Repetition_${repetition}/Traces

for (( level=2 ; level <= 729 ; level++ ))
	do
        echo ~/bin/OptimTiles2D -t ${nbthreads} -n 1024 -i ../Optimize_L2Discrepancy_For_2D_Tiles/Data/Optim_Input/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions/Repetition_${repetition}/Output/2D_0m2net_set_1_level_${level}_Opt.dat
        ~/bin/OptimTiles2D -t ${nbthreads} -n 1024 -i ../Optimize_L2Discrepancy_For_2D_Tiles/Data/Optim_Input/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions/Repetition_${repetition}/Output/2D_0m2net_set_1_level_${level}_Opt.dat | tee ../Repetitions/Repetition_${repetition}/Traces/Trace_Level_${level}.dat
	done
