#!/bin/bash

if [ $# -ne 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

repetition=$1



nthreads=32

for (( level=243; level <= 729; level++ ))
	do
		sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions/Repetition_${repetition}/Traces/slurm-%j-${level}.out" runOneRepetitionOneLevel.sh ${repetition} ${nthreads} ${level}
	done
