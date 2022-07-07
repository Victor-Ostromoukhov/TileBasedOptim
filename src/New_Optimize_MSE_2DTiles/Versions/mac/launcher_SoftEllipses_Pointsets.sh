#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

ind=$1
dir=$2
suffix=$3

nthreads=8

./oneJob_SoftEllipses_Pointsets.sh ${ind} ${nthreads} ${dir} ${suffix}
    
