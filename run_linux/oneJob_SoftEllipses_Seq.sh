#!/bin/bash


ind=$1
nbthreads=$2
fname=$3
suffix=$4

nIterations=128
nItegrandsPerIteration=65536


InputDir="../Tiles_Seq_${suffix}/"
OutputDir="../Output/Tiles_Seq_${suffix}/"
TracesDir="../Traces/Tiles_Seq_${suffix}/"
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}
 
infname=${InputDir}${fname}
outfname=${OutputDir}${fname}
echo cp ${infname} ${outfname}
cp ${infname} ${outfname}

# counters...
#lst=(1 3    9    27    81    243    729   2187   6561 19683 59049)
lst=(1 3    9    27    81    243    729   2187   6561)

lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

for level in $(eval echo {1..$((${lst_length} - 1))})
do
    from=${lst[$((${level} - 1 ))]}
    to=${lst[$level]}
    limit=$((${lst[$((${level} - 1 ))]} + 1))
    infname=${OutputDir}${fname}
    outfname=${OutputDir}${fname}
    echo ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outfname --nbPoints $to --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit
    ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i $infname -o $outfname --outputNextStep $outfname --nbPoints $to --integrandType ${integrandType} -g $nItegrandsPerIteration --limit $limit >> ${TracesDir}/t_${fname}.txt
done



