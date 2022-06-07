#!/bin/bash

if [ $# -ne 2 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads>"
	exit
fi

repetition=$1
nbthreads=$2

mkdir -p Repetition_${repetition}/Output
mkdir -p Repetition_${repetition}/Traces

for (( level=2 ; level <= 242 ; level++ ))
	do
		~/bin/Optimize_L2Discrepancy_For_2D_Tiles -t ${nbthreads} -n 1024 -i ../Optimize_L2Discrepancy_For_2D_Tiles/Data/Optim_Input/2D_0m2net_set_1_level_${level}.dat -o Repetition_${repetition}/Output/2D_0m2net_set_1_level_${level}_Opt.dat >> Repetition_${repetition}/Traces/Trace_Level_${level}.dat
	done
