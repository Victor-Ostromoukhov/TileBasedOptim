#!/bin/bash

nthreads=8

#ntrials=1024 4096 0r 16384 or 65536
ntrials=1024

for (( level=180; level<=199; level++ )) do
    echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials}
    ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} >> traces/tr_${level}.txt
done


