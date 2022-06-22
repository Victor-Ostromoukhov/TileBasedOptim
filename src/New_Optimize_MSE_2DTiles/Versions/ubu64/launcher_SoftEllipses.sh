#!/bin/bash

if [ $# -le 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

repetition=$1
fname=$2


nIterations=256
nItegrandsPerIteration=4096

nthreads=64

lst=(3    9    27    81    243    729)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses


./oneJob_SoftEllipses.sh ${repetition} ${nthreads} ${fname} 2 3 1

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    ./oneJob_SoftEllipses.sh ${repetition} ${nthreads} ${fname} $from $to $limit
done

