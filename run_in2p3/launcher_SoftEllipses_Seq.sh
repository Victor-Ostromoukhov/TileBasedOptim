#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
fname=$2
continueFlag=$3
nthreads=$4

sbatch -n ${nthreads} -J S${ind} -o "slurm-%j-${fname}.out"  ./oneJob_SoftEllipses_Seq.sh ${ind} ${nthreads} ${fname} ${continueFlag}
    
