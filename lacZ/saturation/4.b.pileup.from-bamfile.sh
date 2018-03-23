#!/bin/sh

# Loop through each bam file
for bamfile in ./*.bam
do

# grab filename
fname=${bamfile%.bam}
fname=${bamfile##*/}
echo "Starting screen for $fname"

# Create screen using filename
# screen -d -m -S "$fname".

# Grab date in format dd.mm.yyyy
# now="$(date +'%d.%m.%Y')"

# Run samtools pileup on each
screen -d -m -S "$fname". bash -c "cat $bamfile | samtools mpileup -B -Q 0 -d 500000 -f /home/ngs/references/lacz/lacZFlncCorr/ReferenceGenome.fasta - > $fname'.pileUp.pu' "

done

