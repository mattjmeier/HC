#!/bin/bash

# Loop through chromosomes
echo "What folder would you like to compute coverage for?"
read -e rootdir
echo "Where should the bedgraph files be written?"
read -e outdir

for i in $(find $rootdir -name "*.bam")
do

#for chr in $(cat $mm10_chr)
#for chr in $(echo chrM)
#do
#echo "Starting for $chr on $file"


fname=${i##*/}
fname=${fname%.*}
echo $i
echo $fname

 /utils/appz/samtools/samtools-1.3/samtools view -F 0x400 -b $i | genomeCoverageBed -bga -ibam stdin -g /s0/ngs/references/mouse/ucsc/mm10/index/samtools/mm10.genome > $outdir/$fname.bga &

# Don't start more than X instances
#numProc=$(ps auxw|grep -i genomeCoverageBed|grep -v grep|wc -l)
#while [ $numProc -gt 15 ]
#do sleep 5
#numProc=$(ps auxw|grep -i genomeCoverageBed|grep -v grep|wc -l)
#done

done


