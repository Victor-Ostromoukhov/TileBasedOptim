#!/bin/bash

nthreads=$1
level=$2

srcdir=optim_data_2D
resdir=optim_output_2D

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

#ntrials=1024 4096 0r 16384 or 65536
ntrials=16384

echo ~/bin/OptimTiles2D -t ${nthreads} -i ${srcdir}/2D_0m2net_set_1_level_${level}.dat -o ${resdir}/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} > traces/tr_${level}.txt
~/bin/OptimTiles2D -t ${nthreads} -i ${srcdir}/2D_0m2net_set_1_level_${level}.dat -o ${resdir}/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} >> traces/tr_${level}.txt
