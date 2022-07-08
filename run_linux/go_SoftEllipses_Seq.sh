#!/bin/bash

if [ $# -lt 3 ] ;
then
    echo "Usage : go_SoftEllipses_Seq <set_from> <set_to> <nthreads>"
    exit
fi

set_from=$1
set_to=$2
nthreads=$3

suffix="PrevLevel"

inputFiles=(`ls ../Tiles_Seq_${suffix}/`)

lst_length=${#inputFiles[@]}
#lst_length=1


for (( ind=set_from ; ind <= set_to ; ind++ ))
	do
        fname=${inputFiles[$((${ind} - 1))]}
        ./launcher_SoftEllipses_Seq.sh $ind ${fname} ${suffix} ${nthreads}
	done
