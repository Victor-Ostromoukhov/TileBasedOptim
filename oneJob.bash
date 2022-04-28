#!/bin/bash

nthreads=$1
level=$2

echo ~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat > traces/tr_${level}.txt
~/bin/OptimTiles2D -t ${nthreads} -i optim_data/2D_0m2net_set_1_level_${level}.dat -o optim_output/2D_0m2net_set_1_level_${level}.dat >> traces/tr_${level}.txt
