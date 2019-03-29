#!/usr/bin/env bash
cpu=$(lscpu | awk 'NR==4' | awk '{print $2}')
#number of cpus

bucketDir=s3://deepmod-benchmark/miten-test-data/
#bucket

fastq=($(aws s3 ls ${bucketDir} | grep "fastq" | awk '{print $4}'))
fast5=($(aws s3 ls ${bucketDir} | grep "fast5" | awk '{print $4}'))
#sequencing files

for i in $( eval echo {0..$(($cpu-1))} ); do touch S$i-fastq.txt; touch S$i-fast5.txt; done
#create file tables

for i in $( eval echo {0..$((${#fastq[@]}-1))..$cpu} ); do
for j in $( eval echo {0..$(($cpu-1))} ); do
echo ${fastq[(($i+$j))]} >> S$j-fastq.txt
echo ${fast5[(($i+$j))]} >> S$j-fast5.txt
done
done
#distribute reads

for i in $( eval echo {0..$(($cpu-1))} ); do screen -S S$i; done 
#create screens

sh workflow.sh
#initiate workflow in each screen

pid=($(ps -A | grep "screen" | awk '{print $1}'))
#pid of the screens

for i in $( eval echo {0..$(($cpu-1))} ); do taskset -cp $i ${pid[i]}; done
#parallel screens

top
1
#monitoring cpus

pkill screen
#kill all screens
