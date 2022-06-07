#!/bin/bash

if [ $# -ne 3 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

repetition=$1
nbthreads=$2
level=$3

mkdir -p ../Repetitions/Repetition_${repetition}/Output
mkdir -p ../Repetitions/Repetition_${repetition}/Traces

ntrials=1024

../Optimize_L2Discrepancy_For_2D_Tiles/build/bin/Optimize_L2Discrepancy_For_2D_Tiles -t ${nbthreads} -n ${ntrials} -i ../Optimize_L2Discrepancy_For_2D_Tiles/Data/Optim_Input/2D_0m2net_set_1_level_${level}.dat -o ../Repetitions/Repetition_${repetition}/Output/2D_0m2net_set_1_level_${level}_Opt.dat
