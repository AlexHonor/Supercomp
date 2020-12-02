#!/bin/bash -x
module load SpectrumMPI
mpixlC -O3 -lm -std=c++98 task2.cpp -o task2

rm -f PolusResults.csv
printf "ProcNum,TotalTime,IterationNumber,CellsNumber,HighBorder,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,TimeSpent\n" > PolusResults.csv

for procNum in 1 10 20 40
do
    for cellSize in 128 256 512
    do 
        for highBorder in 1.0 3.1415
        do 
            mpisubmit.pl -p $procNum ../task2 -- $cellSize 0.025 20 $highBorder PolusResults.csv
        done
    done
done