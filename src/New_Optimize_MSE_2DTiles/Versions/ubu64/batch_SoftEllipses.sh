#!/bin/bash


inputFiles=(`ls ../../Data/Input/Tiles/`)

lst_length=${#inputFiles[@]}


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        fname=${inputFiles[$((${ind}))]}
        ./runOneRepetitionLauncher_SoftEllipses.sh $ind $fname
	done
