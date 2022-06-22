#!/bin/bash

for (( i = 2; i < 729; i++ )); do
      mv 2D_0m2net_set_1_level_${i}.dat 3D_0m2net_set_1_level_${i}.dat
done
