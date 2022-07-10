#!/bin/bash
# counters...
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)

suffix="PrevLevel"

ind=$1
nbthreads=$2
fname=$3
continueFlag=$4

echo continueFlag=$continueFlag


if [ $continueFlag = false ]; then
    InputDir="../Tiles_Seq_${suffix}/"
    lst=(1 3    9    27    81    243    729   2187   6561)
else
    InputDir="../Output/Tiles_Seq_${suffix}/"
    lst=(729   2187   6561)
fi
OutputDir="../Output/Tiles_Seq_${suffix}/"
TracesDir="../Traces"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
infname=${InputDir}${fname}
outfname=${OutputDir}${fname}
if [ $continueFlag = false ]; then
    echo cp ${infname} ${outfname}
    cp ${infname} ${outfname}
fi

echo InputDir=$InputDir
echo OutputDir=$OutputDir
echo TracesDir=$TracesDir

nIterations=128
nItegrandsPerIteration=65536

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    npts=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    infname=${OutputDir}${fname}
    outfname=${OutputDir}${fname}
    if [ ${npts} -gt 729 ]; then
        nIterations=$((${nIterations} * 2 ))
    fi
    echo trace into ${TracesDir}/t_${fname}.txt
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints $npts --limit $limit -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration 
    echo ~/bin/Optimize_MSE_2DTiles --nbPoints $npts --limit $limit -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration >> ${TracesDir}/t_${fname}.txt
    ~/bin/Optimize_MSE_2DTiles --nbPoints $npts --limit $limit -t ${nbthreads} -n $nIterations -i $infname -o $outfname --integrandType ${integrandType} -g $nItegrandsPerIteration >> ${TracesDir}/t_${fname}.txt
done



