#!/bin/bash


inputFiles=(`ls ../../Data/Input/Tiles_Seq/`)

lst_length=${#inputFiles[@]}

for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        fname=${inputFiles[$((${ind} - 1))]}
       ./launcher_SoftEllipses_Seq.sh $ind $fname
	done
