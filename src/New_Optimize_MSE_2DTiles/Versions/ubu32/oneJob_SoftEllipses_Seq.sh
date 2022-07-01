#!/bin/bash

if [ $# -le 2 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads> <level>"
	exit
fi

ind=$1
nbthreads=$2
fname=$3

nIterations=16
nItegrandsPerIteration=65536


InputDir="../../Data/Input/Tiles_Seq_PrevLevel/"
OutputDir="../../Data/Output/"
TracesDir="../../Data/Traces"
mkdir -p ${NextIterDir}
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
 
    infname=${InputDir}${fname}
    outfname=${OutputDir}${fname}
echo cp ${infname} ${outfname}
cp ${infname} ${outfname}

#the rest of counters...
#lst=(1 3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
#lst=(1 3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)
lst=(6    9    13    19    27    39    56    81    117    168    243    350    505    729)
#lst=(1 3    9    27    81    243    729)
#lst=(27 81 243)

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    infname=${OutputDir}${fname}
    outfname=${OutputDir}${fname}
    ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outfname --nbPoints $to --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit
done



