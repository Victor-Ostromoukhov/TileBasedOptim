#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
fname=$2
suffix=$3

nthreads=32


#echo "Starting Job" > ../../Data/Traces/trace-${fname}.txt
sbatch -n ${nthreads} -J S${ind} -o "../../Data/Traces/slurm-%j-${fname}.out" oneJob_SoftEllipses_Seq.sh ${ind} ${nthreads} ${fname} ${suffix}
    
