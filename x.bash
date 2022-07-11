#!/bin/bash

if [ `uname` = "Darwin" ]; then
    mexec="/Applications/Mathematica.app/Contents/MacOS/MathKernel "
else
    mexec="/usr/local/bin/math"
fi


if [ -n "$1" ]
then
    nproc=$1
else
    nproc=1
fi

cp TileBasedOptim/TileBasedOptim.m tmp/proc
echo " " >> tmp/proc
echo "prepSoftEllipses2D[9,9]" >> tmp/proc
chmod 700 tmp/proc
$mexec < tmp/proc > tmp/trace_prepSoftEllipses2D.txt &

echo prepSoftEllipses2D[9,9] has been launched
echo watch tmp/trace_prepSoftEllipses2D.txt
