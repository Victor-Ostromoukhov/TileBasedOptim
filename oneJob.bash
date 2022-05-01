#!/bin/bash

nthreads=$1
level=$2

#ntrials=4096 0r 16384 or 65536
ntrials=16384

echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} > traces/tr_${level}.txt
~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} >> traces/tr_${level}.txt
