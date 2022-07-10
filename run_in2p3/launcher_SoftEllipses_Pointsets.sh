#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
dir=$2
continueFlag=$3
nthreads=$4

sbatch -n ${nthreads} -J P${ind} -o "slurm-%j-${dir}.txt" ./oneJob_SoftEllipses_Pointsets.sh ${ind} ${nthreads} ${dir} ${continueFlag}
    
