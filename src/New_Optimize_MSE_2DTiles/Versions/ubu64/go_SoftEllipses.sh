#!/bin/bash


inputFiles=(`ls ../../Data/Input/Tiles/`)

lst_length=${#inputFiles[@]}
lst_length=1


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        fname=${inputFiles[$((${ind}))]}
       ./launcher_SoftEllipses.sh $ind $fname
	done
