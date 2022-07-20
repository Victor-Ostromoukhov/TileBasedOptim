#!/bin/bash
suffix="VarSize"

#if [ $# -lt 1 ] ;
#then
#    echo "Usage : go_SoftEllipses_Pointsets <set_from> <set_to> <nthreads> [<continueFlag>]"
#    exit
#fi
#
#set_from=$1
#set_to=$2
#nthreads=$3
#if [ $# -ge 4 ]; then
#    continueFlag=$4
#else
#    continueFlag=false
#fi
#

nthreads=64
continueFlag=false

inputDirs=(`ls ../Tiles_Pointsets_${suffix}/`)

lst_length=${#inputDirs[@]}
#lst_length=1


for (( ind=1 ; ind <= lst_length ; ind++ ))
	do
        dir=${inputDirs[$((${ind} - 1))]}
       ./launcher_SoftEllipses_Pointsets.sh $ind ${dir} ${continueFlag}  ${nthreads}
	done
