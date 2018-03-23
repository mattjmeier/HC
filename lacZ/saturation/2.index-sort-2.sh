#!/bin/sh

# Loop through each bam file
for samfile in ./*.sam

do

# grab filename
fname=${samfile%.sam}
fname=${fname##*/}
echo "Starting screen for $fname"

# Create screen using filename
# screen -d -m -S "$fname".

# Grab date in format dd.mm.yyyy
# now="$(date +'%d.%m.%Y')"

# Run samtools view and sort and index on each
screen -d -m -S "$fname". bash -c "samtools view -bh $samfile | samtools sort -m 35G - > $fname.bam"

# If the number of samtools processes that are running exceeds a certain number, sleep for 5 seconds, then check again until a process ends.

numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)

while [ $numProc -gt 8 ]
do sleep 5
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)

done

done
