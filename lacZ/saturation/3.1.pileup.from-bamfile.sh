#!/bin/sh

# Loop through each bam file
for bamfile in $(echo "./*_015*.bam" "./*_036*.bam")
do

# grab filename
fname=${bamfile%.bam}
fname=${fname##*/}
echo "Starting screen for $fname"

# Create screen using filename
screen -d -m -S "$fname".

# Grab date in format dd.mm.yyyy
# now="$(date +'%d.%m.%Y')"

# Run samtools pileup on each
screen -S "$fname". -p 0 -X screen bash -c " samtools mpileup -B -Q 0 -d 500000 -f /s0/ngs/references/lacz/lacZFlncCorr/ReferenceGenome.fasta $bamfile > $fname'.pileUp.pu' "


# If the number of bowtie processes that are running exceeds a certain number, sleep for 5 seconds, then check again until a process ends.

numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)

while [ $numProc -gt 8 ]
do sleep 1
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)

done

done

