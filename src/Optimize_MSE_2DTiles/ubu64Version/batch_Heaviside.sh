#!/bin/bash

for (( rep=1 ; rep <= 64 ; rep++ ))
	do
        ./runOneRepetition_Heaviside.sh $rep 64
	done
