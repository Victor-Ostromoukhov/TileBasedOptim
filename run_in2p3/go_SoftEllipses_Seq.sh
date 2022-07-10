#!/bin/bash
suffix="PrevLevel"

if [ $# -lt 3 ] ;
then
    echo "Usage : go_SoftEllipses_Seq <set_from> <set_to> <nthreads> [<continueFlag>]"
    exit
fi

set_from=$1
set_to=$2
nthreads=$3
if [ $# -ge 4 ]; then
    continueFlag=$4
else
    continueFlag=false
fi

inputFiles=(`ls ../Tiles_Seq_${suffix}/`)

lst_length=${#inputFiles[@]}
#lst_length=1


for (( ind=set_from ; ind <= set_to ; ind++ ))
	do
        fname=${inputFiles[$((${ind} - 1))]}
        ./launcher_SoftEllipses_Seq.sh $ind ${fname} ${continueFlag} ${nthreads}
	done
