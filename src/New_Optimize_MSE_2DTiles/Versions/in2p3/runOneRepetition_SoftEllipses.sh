
#!/bin/bash

if [ $# -ne 2 ] ;
then
	echo "Usage : runOneJob <repetition> <nbthreads>"
	exit
fi

if [ $# -ge 1 ]; then
  repetition=$1
else
  repetition=99
fi

if [ $# -ge 2 ]; then
  nbthreads=$2
else
  nbthreads=64
fi
InputDir="../../Data/Input/Tiles/"
NextIterDir="../../Data/Input/TilesLimited/"
OutputDir="../../Data/Output/Limited/SoftEllipses/OutputTiles/Repetition_${repetition}/"
TracesDir="../../Data/Output/Limited/SoftEllipses/Traces/Repetition_${repetition}/"
nIterations=256
#nItegrandsPerIteration=16384
nItegrandsPerIteration=4096
#lst=(3 4 5 6 7 8 9 10 11 13 15 17 19 21 24 27 31 34 39 44 50 56 63 72 81 92 103 117 132 149 168 190 215 243 275 310 350 396 447 505 571 645 729)
#lst=(3    4    6    9    13    19    27    39    56    81    117    168    243    350    505    729)

lst=(3    9    27    81    243    729)
lst_length=${#lst[@]}
integrandType=2 # SoftEllipses

mkdir -p ${NextIterDir}
mkdir -p ${OutputDir}
mkdir -p ${TracesDir}

  ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n ${nIterations} -i "${InputDir}2D_0m2net_set_1_level_729.dat" -o "${OutputDir}2D_0m2net_set_1_level_Opt${lst[0]}.dat" --outputNextStep "${NextIterDir}2D_0m2net_set_1_level_OptNext${lst[0]}.dat" --nbPoints ${lst[0]} --integrandType ${integrandType} -g ${nItegrandsPerIteration} -l 1



for level in $(eval echo {1..$((${lst_length} - 1))})
  do
#echo        ~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i "${NextIterDir}2D_0m2net_set_1_level_OptNext$((${lst[$((${level} - 1 ))]})).dat" -o "${OutputDir}2D_0m2net_set_1_level_Opt${lst[$level]}.dat" --outputNextStep "${NextIterDir}2D_0m2net_set_1_level_OptNext$((${lst[$((${level} - 1 ))]} + 1)).dat" --nbPoints ${lst[$level]} --integrandType ${integrandType} -g $nItegrandsPerIteration -l $((${lst[$((${level} - 1 ))]} + 1)) | tee ../Repetitions_Heaviside_ubu32/Repetition_${repetition}/Limited/Traces/Trace_Level_${level}.dat
~/bin/Optimize_MSE_2DTiles -t ${nbthreads} -n $nIterations -i "${NextIterDir}2D_0m2net_set_1_level_OptNext$((${lst[$((${level} - 1 ))]})).dat" -o "${OutputDir}2D_0m2net_set_1_level_Opt${lst[$level]}.dat" --outputNextStep "${NextIterDir}2D_0m2net_set_1_level_OptNext$((${lst[$((${level}))]})).dat" --nbPoints ${lst[$level]} --integrandType ${integrandType} -g $nItegrandsPerIteration -l $((${lst[$((${level} - 1 ))]} + 1)) | tee "${TracesDir}Trace_Level_${level}.dat"
	done

