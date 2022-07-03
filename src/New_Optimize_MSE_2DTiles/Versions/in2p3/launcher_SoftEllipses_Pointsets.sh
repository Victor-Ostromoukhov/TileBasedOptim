#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
dir=$2
suffix=$3

nthreads=32

sbatch -n ${nthreads} -J R${ind} -o "../../Data/Traces/slurm-%j-${fname}.out" ./oneJob_SoftEllipses_Pointsets.sh ${ind} ${nthreads} ${dir} ${suffix}
    
