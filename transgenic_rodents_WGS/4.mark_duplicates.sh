#!/bin/bash

files=$(find . -maxdepth 1 -name "*sort.merge.bam")
for i in $files
do

#java -jar /utils/appz/picard/picard-tools-1.96/SortSam.jar INPUT=$i OUTPUT=$i.sorted_reads.bam SORT_ORDER=coordinate

java -jar /utils/appz/picard/picard/dist/picard.jar MarkDuplicates INPUT=$i OUTPUT=$i.dedup_reads.bam  METRICS_FILE=$i.metrics.txt 

done
