#!/bin/bash

nthreads=32

for (( level=2; s<=10; s++ )) do
   for (( trial=1; trial<=1; trial++ )) do
        sbatch -c ${nthreads} -J L${level} oneJob.bash ${nthreads} ${level}
    done
done


