#!/bin/bash

nthreads=64

#ntrials=16384 or 65536
ntrials=16384

for (( level=137; level<=199; level++ )) do
    echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials}
    ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n ${ntrials} >> traces/tr_${level}.txt
done


