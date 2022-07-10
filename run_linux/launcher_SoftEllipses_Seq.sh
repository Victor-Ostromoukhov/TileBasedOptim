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

./oneJob_SoftEllipses_Seq.sh ${ind} ${nthreads} ${fname} ${continueFlag}
    
