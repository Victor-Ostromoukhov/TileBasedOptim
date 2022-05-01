#!/bin/bash

nthreads=32

for (( level=712; level<=712; level++ )) do
         #710
         #sbatch -J L${level} oneJob.bash ${nthreads} ${level}
          #711
         #sbatch -c ${nthreads} -J L${level} oneJob.bash ${nthreads} ${level}
          #712
         sbatch -n ${nthreads} -J L${level} oneJob.bash ${nthreads} ${level}
done


