#!/bin/bash

if [ $# -lt 3 ] ;
then
    echo "Usage : go_SoftEllipses_Pointsets <set_from> <set_to> <nthreads>"
    exit
fi

set_from=$1
set_to=$2
nthreads=$3

suffix="PrevLevel"

inputDirs=(`ls ../Tiles_Pointsets_${suffix}/`)

lst_length=${#inputDirs[@]}
#lst_length=1


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        dir=${inputDirs[$((${ind} - 1))]}
       ./launcher_SoftEllipses_Pointsets.sh $ind ${dir} ${suffix}  ${nthreads}
	done
