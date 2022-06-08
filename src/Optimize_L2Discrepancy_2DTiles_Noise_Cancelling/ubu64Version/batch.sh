#!/bin/bash

for (( rep=6 ; rep <= 64 ; rep++ ))
	do
        ./runOneRepetition.sh $rep 64
	done
