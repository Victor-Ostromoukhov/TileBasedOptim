#!/bin/bash

if [ $# -lt 1 ] ;
then
    echo "Usage : go_SoftEllipses_Seq <suffix>            PrevLevel or CurLevel "
    exit
fi

suffix=$1

inputFiles=(`ls ../../Data/Input/Tiles_Seq_${suffix}/`)

lst_length=${#inputFiles[@]}
#lst_length=1


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        fname=${inputFiles[$((${ind} - 1))]}
       ./launcher_SoftEllipses_Seq.sh $ind ${fname} ${suffix}
	done
