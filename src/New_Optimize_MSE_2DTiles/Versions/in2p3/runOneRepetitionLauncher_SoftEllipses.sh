#!/bin/bash

if [ $# -ne 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

repetition=$1


nIterations=256
#nItegrandsPerIteration=16384
nItegrandsPerIteration=4096

nthreads=32


lst=(3    9    27    81    243    729)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses


echo sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/slurm-%j-${level}.out" OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads} 2 3 1
OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads} 2 3 1

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    echo sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/slurm-%j-${level}.out" OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads}  $from $to $limit
   #./OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads}  $from $to $limit
done

