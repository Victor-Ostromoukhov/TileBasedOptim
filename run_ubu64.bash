#!/bin/bash

nthreads=64

for (( level=2; level<=729; level++ )) do
        echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n 100000 > traces/tr_${level}.txt
        ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat -n 100000 >> traces/tr_${level}.txt
done


