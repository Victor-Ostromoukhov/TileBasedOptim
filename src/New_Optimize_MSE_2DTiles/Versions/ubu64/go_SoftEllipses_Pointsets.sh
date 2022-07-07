#!/bin/bash

#if [ $# -lt 1 ] ;
#then
#    echo "Usage : go_SoftEllipses_Seq <suffix>            PrevLevel or CurLevel "
#    exit
#fi

suffix="PrevLevel"

inputDirs=(`ls ../../Data/Input/Tiles_Pointsets_${suffix}/`)

lst_length=${#inputDirs[@]}
#lst_length=1


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        dir=${inputDirs[$((${ind} - 1))]}
       ./launcher_SoftEllipses_Pointsets.sh $ind ${dir} ${suffix}
	done
