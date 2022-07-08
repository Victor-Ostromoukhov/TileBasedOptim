#!/bin/bash


ind=$1
nbthreads=$2
dir=$3
suffix=$4

nIterations=128
nItegrandsPerIteration=65536

InputDir="../Tiles_Pointsets_${suffix}/${dir}/"
OutputDir="../Output/Tiles_Pointsets_${suffix}/"
TracesDir="../Traces/Tiles_Pointsets_${suffix}/"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
 
# counters...
#lst=(1 3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
#lst=(1 3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)
#lst=(3    4    5    6    7    8    9    10    11    13    15    17    19    21    24    27    31    34    39    44    50    56    63    72    81    92    103    117    132    149    168    190    215    243    275    310    350    396    447    505    571    645    729    824    931    1051    1188    1342    1516    1713    1936    2187    2471    2792    3154    3564    4026    4549    5140    5807    6561)
#lst=(1 3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)
lst=(1 3    9    27    81    243    729)

inputFiles=(`ls ${InputDir}`)

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for (( ind=0 ; ind < lst_length ; ind++ ))
do
    npts=${lst[$((${ind}+1))]}
    infname=${InputDir}${inputFiles[$((${ind}))]}
    outfname=${OutputDir}${inputFiles[$((${ind}))]}
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration
    ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration  | tee ${TracesDir}/t_${fname}.txt
done



