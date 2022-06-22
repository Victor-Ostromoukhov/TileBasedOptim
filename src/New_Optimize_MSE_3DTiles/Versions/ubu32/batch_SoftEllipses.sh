#!/bin/bash

for (( rep=3 ; rep <= 16 ; rep++ ))
	do
        ./runOneRepetition_SoftEllipses.sh $rep 32
	done
