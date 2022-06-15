#!/bin/bash

for (( rep=2 ; rep <= 16 ; rep++ ))
	do
        ./runOneRepetition_SoftEllipses.sh $rep 64
	done
