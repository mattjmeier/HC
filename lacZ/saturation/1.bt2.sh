#!/bin/sh

# Loop through each bam file

for fqfile in ./*.fastq
do

# grab filename
fname=${fqfile%.fastq}
fname=${fname##*/}

# Create screen using filename
echo "Starting screen for $fname"
screen -d -m -S "$fname".

# Run bt2 on each
screen -S "$fname". -p 0 -X screen bash -c "/utils/appz/bowtie2/bowtie2-2.1.0/bowtie2 --very-sensitive-local -x /s0/ngs/references/lacz/lacZFlncCorr/lacZFlncCorr -p 4 -q $fqfile --rg-id id --rg PL:iontorrent -S $fname.sam &> $fname.stats"

# If the number of bowtie processes that are running exceeds a certain number, sleep for 5 seconds, then check again until a process ends.

numProc=$(ps auxw|grep -i bowtie|grep -v grep|wc -l)

while [ $numProc -gt 7 ]
do sleep 5
numProc=$(ps auxw|grep -i bowtie|grep -v grep|wc -l)

done


done

