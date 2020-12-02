#!/bin/bash -x

mpicxx task2.cpp -o task2

mpisubmit.pl -p 10 -w 00:15 task2 -- test.ini

#for procNum in 10 20 40
#do 
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_128.ini
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_128p.ini
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_256.ini
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_256p.ini
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_512.ini
#	mpisubmit.pl -p $procNum -w 00:15 main -- configs/config_512p.ini
#done