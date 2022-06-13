#!/bin/bash

for (( rep=3 ; rep <= 16 ; rep++ ))
	do
        ./runOneRepetition_Heaviside.sh $rep 32
	done
