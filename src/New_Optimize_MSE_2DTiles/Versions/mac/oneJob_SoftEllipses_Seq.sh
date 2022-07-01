#!/bin/bash

if [ $# -le 2 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

ind=$1
nbthreads=$2
fname=$3

nIterations=8

#nItegrandsPerIteration=4096
nItegrandsPerIteration=16384
#nItegrandsPerIteration=262144


InputDir="../../Data/Input/Tiles_Seq/"
NextIterDir="../../Data/Output/"
OutputDir="../../Data/Output/"
TracesDir="../../Data/Traces"
mkdir -p ${NextIterDir}
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
  
#first run for 3 pts
from=2
to=3
limit=1
    infname=${InputDir}${fname}
    outfname=${OutputDir}${fname}
    outputNextStep=${NextIterDir}${fname}

~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outputNextStep --nbPoints $to --integrandType 2 -g $nItegrandsPerIteration --limit $limit

#the rest of counters...
#lst=(3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
lst=(3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)
#lst=(3    9    27    81    243    729)

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    infname=${NextIterDir}${fname}
    outfname=${OutputDir}${fname}
    outputNextStep=${NextIterDir}${fname}
    ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outputNextStep --nbPoints $to --integrandType 2 -g $nItegrandsPerIteration --limit $limit
done



