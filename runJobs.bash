#!/bin/bash

nthreads=32

for (( level=2; level<=729; level++ )) do
   for (( trial=1; trial<=1; trial++ )) do
        sbatch -c ${nthreads} -J L${level} oneJob.bash ${nthreads} ${level}
    done
done


