#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
fname=$2

nIterations=256
nItegrandsPerIteration=4096

nthreads=8

./oneJob_SoftEllipses_Seq.sh ${ind} ${nthreads} ${fname}
    
