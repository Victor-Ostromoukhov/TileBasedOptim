#!/bin/bash
suffix="PrevLevel"

ind=$1
nbthreads=$2
dir=$3
continueFlag=$4

if [ $continueFlag ] ;
then
    InputDir="../Output/Tiles_Pointsets_${suffix}/${dir}/"
else
    InputDir="../Tiles_Pointsets_${suffix}/${dir}/"
fi
OutputDir="../Output/Tiles_Pointsets_${suffix}/${dir}/"
TracesDir="../Traces/Tiles_Pointsets_${suffix}/${dir}/"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}

# counters...
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)
#lst=(1 3    9    27    81    243    729)
#lst=(1 3    9    27    81    243    729   2187   6561)
lst=(729   2187   6561)

nIterations=128
nItegrandsPerIteration=65536

inputFiles=(`ls ${InputDir}`)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for (( ind=0 ; ind < lst_length ; ind++ ))
do
    npts=${lst[$((${ind}+1))]}
    fname=${inputFiles[$((${ind}))]}
    infname=${InputDir}${fname}
    outfname=${OutputDir}${fname}
    if [ ${npts} -gt 729 ] ;
    then
        nIterations=$((${nIterations} * 2 ))
    fi
echo ${npts}
echo ${fname}
echo ${infname}
echo ${outfname}
echo ${nIterations}
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration  >> ${TracesDir}/t_${fname}.txt
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints ${npts} -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration  >> ${TracesDir}/t_${fname}.txt
done


