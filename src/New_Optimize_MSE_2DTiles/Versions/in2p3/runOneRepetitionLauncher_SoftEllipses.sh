#!/bin/bash

if [ $# -ne 1 ] ;
then
	echo " Usage : runJobs <repetition>"
	exit
fi

repetition=$1

InputDir="../../Data/Input/Tiles/"
NextIterDir="../../Data/Input/TilesLimited/"
OutputDir="../../Data/Output/Limited/SoftEllipses/OutputTiles/Repetition_${repetition}/"
TracesDir="../../Data/Output/Limited/SoftEllipses/Traces/Repetition_${repetition}/"
nIterations=256
#nItegrandsPerIteration=16384
nItegrandsPerIteration=4096

nthreads=32

###for (( level=243; level <= 729; level++ ))
#lst=(3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
#lst=(3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)
#for level in ${lst[@]}
#	do
#		sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/slurm-%j-${level}.out" OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads} ${level}
#	done


lst=(3    9    27    81    243    729)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

mkdir -p ${NextIterDir}
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
  
sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/slurm-%j-${level}.out" OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads} 2 3 1

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    sbatch -n ${nthreads} -J R${repetition}L${level} -o "../Repetitions_SoftEllipses/Repetition_${repetition}/Traces/slurm-%j-${level}.out" OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads}  $from $to $limit
   ./OneRepetitionOneLevel_SoftEllipses.sh ${repetition} ${nthreads}  $from $to $limit
done

