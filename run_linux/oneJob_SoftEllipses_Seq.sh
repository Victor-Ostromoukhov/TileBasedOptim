#!/bin/bash
# counters...
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)

suffix="PrevLevel"

ind=$1
nbthreads=$2
fname=$3
continueFlag=$4

if [ continueFlag ] ;
then
    InputDir="../Output/Tiles_Seq_${suffix}/"
    lst=(729   2187   6561)
else
    InputDir="../Tiles_Seq_${suffix}/"
    lst=(1 3    9    27    81    243    729   2187   6561)
fi
OutputDir="../Output/Tiles_Seq_${suffix}/"
TracesDir="../Traces/Tiles_Seq_${suffix}/"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
 
infname=${InputDir}${fname}
outfname=${OutputDir}${fname}
echo cp ${infname} ${outfname}
cp ${infname} ${outfname}

nIterations=64
nItegrandsPerIteration=65536

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    if [ level -gt 729 ] ;
    then
        nIterations=$((${nIterations} * 2 ))
    fi
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    infname=${OutputDir}${fname}
    outfname=${OutputDir}${fname}
    echo ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outfname --nbPoints $to --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit
    ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outfname --nbPoints $to --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit >> ${TracesDir}/t_${fname}.txt
done



