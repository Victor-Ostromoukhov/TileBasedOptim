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

sbatch -n ${nthreads} -J R${ind} -o "../../Traces/slurm-%j-${fname}.out" oneJob_SoftEllipses_Seq.sh ${ind} ${nthreads} ${fname}  ${suffix}
    
