#!/bin/bash

nthreads=8

for (( level=200; level<=729; level++ )) do
    echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n 1024 
    ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n 1024 >> traces/tr_${level}.txt
done


