#!/bin/bash

if [ $# -ne 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

repetition=$1

mkdir -p ../Repetitions_Heaviside/Repetition_${repetition}/Output
mkdir -p ../Repetitions_Heaviside/Repetition_${repetition}/Traces

nthreads=32

###for (( level=243; level <= 729; level++ ))
lst=(3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
for level in ${lst[@]}
	do
		sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_Heaviside/Repetition_${repetition}/Traces/slurm-%j-${level}.out" runOneRepetitionOneLevel.sh ${repetition} ${nthreads} ${level}
	done
