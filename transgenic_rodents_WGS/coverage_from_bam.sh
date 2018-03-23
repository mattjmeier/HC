#!/bin/bash

# Loop through chromosomes
echo "What file would you like to compute coverage for?"
read -e file
echo "Where should the coverage files be written?"
read -e outdir

for chr in $(cat $mm10_chr)
#for chr in $(echo chrM)
do
echo "Starting for $chr on $file"

fname=${file##*/}
fname=${fname%.*}
echo $file
echo $fname

/utils/appz/samtools/samtools-1.3/samtools view -F 0x400 -b $file $chr | genomeCoverageBed -ibam stdin -g /s0/ngs/references/mouse/ucsc/mm10/index/samtools/mm10.genome > $outdir/$fname.$chr.COVERAGE_HISTOGRAM.bed &

# Don't start more than X instances
numProc=$(ps auxw|grep -i genomeCoverageBed|grep -v grep|wc -l)
while [ $numProc -gt 13 ]
do sleep 5
numProc=$(ps auxw|grep -i genomeCoverageBed|grep -v grep|wc -l)
done



done


