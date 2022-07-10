#!/bin/bash
suffix="PrevLevel"

ind=$1
nbthreads=$2
dir=$3
if [ $# -ge 4 ]; then
     continueFlag=$4
else
    continueFlag=false
fi

echo $continueFlag

if [ $continueFlag ]; then
    InputDir="../Output/Tiles_Pointsets_${suffix}/${dir}/"
    echo $continueFlag
else
    InputDir="../Tiles_Pointsets_${suffix}/${dir}/"
fi
OutputDir="../Output/Tiles_Pointsets_${suffix}/${dir}/"
TracesDir="../Traces/Tiles_Pointsets_${suffix}/${dir}/"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}

echo InputDir=$InputDir
echo OutputDir=$OutputDir
# counters...
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)
#lst=(1 3    9    27    81    243    729)
lst=(1 3    9    27    81    243    729   2187   6561)

nIterations=128
nItegrandsPerIteration=65536

inputFiles=(`ls ${InputDir}`)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for (( ind=6 ; ind < lst_length-1 ; ind++ ))
do
    npts=${lst[$((${ind}+1))]}
    fname=${inputFiles[$((${ind}))]}
    infname=${InputDir}${fname}
    outfname=${OutputDir}${fname}
    if [ ${npts} -gt 729 ] ;
    then
        nIterations=$((${nIterations} * 2 ))
        echo ========= ${nIterations}
    fi
    echo $infname
    echo $outfname
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration  >> ${TracesDir}/t_${fname}.txt
    ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration  >> ${TracesDir}/t_${fname}.txt
done
