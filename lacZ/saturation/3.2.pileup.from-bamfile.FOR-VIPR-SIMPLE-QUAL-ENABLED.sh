#!/bin/sh

# Loop through each bam file
for bamfile in ./*.bam
do

# grab filename
fname=${bamfile%.bam}
fname=${bamfile##*/}
echo "Starting screen for $fname"

# Create screen using filename
screen -d -m -S "$fname".

# Grab date in format dd.mm.yyyy
# now="$(date +'%d.%m.%Y')"

# Run samtools pileup on each
screen -S "$fname". -p 0 -X screen bash -c "cat $bamfile | samtools mpileup -d 1000000 -s -f /home/ngs/references/lacz/lacZFlncCorr/ReferenceGenome.fasta - > $fname'.pileUp.pu'; exit; "


# If the number of samtools processes that are running exceeds a certain number (start with 14), sleep for 5 seconds, then check again until a process ends.
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)
while [ $numProc -gt 35 ]
do sleep 5
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)
done



done

