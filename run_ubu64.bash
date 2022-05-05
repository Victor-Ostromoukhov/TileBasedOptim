#!/bin/bash

nthreads=16

#ntrials=1024 or 4096 or 16384 or 65536 or 262144
ntrials=1024

srcdir=optim_data_2D_MatBuilderOnly
resdir=optim_output_2D_MatBuilderOnly

if [[ ! -d ${srcdir} ]]; then
    echo ${resdir} does not exist
    exit
fi
if [[ ! -d ${resdir} ]]; then
    mkdir ${resdir}
fi
if [[ ! -d traces ]]; then
    mkdir traces
fi

for (( level=2; level<=729; level++ )) do
    echo ~/bin/OptimTiles2D -t ${nthreads} -i ${srcdir}/2D_0m2net_set_1_level_${level}.dat -o ${resdir}/2D_0m2net_set_1_level_${level}.dat -n ${ntrials}
         ~/bin/OptimTiles2D -t ${nthreads} -i ${srcdir}/2D_0m2net_set_1_level_${level}.dat -o ${resdir}/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} >> traces/tr_${level}.txt
done


