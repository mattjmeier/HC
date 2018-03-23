#!/bin/bash
# Takes a directory and samples as input, indexes the BCF files.
# Part one of script - see part 2 to merge.

echo "What directory (or directories, separated by spaces) are your BCF files, split by chromosome, located in?"
read -e rootdir

echo "What sample names should be processed? Please provide a list of exact names, separated by the pipe symbol, |"
read samples

echo "Where should the merged VCF files be written?"
read -e outdir

chromosomes=$(find $rootdir -name "*.bcf" | sed 's/^.*\(chr[0-9A-Za-z]*\)\.bcf/\1/g' | sort | uniq)

for i in $chromosomes
do
echo $i

files=$(find . -regextype posix-egrep -regex ".*($samples).$i.bcf")

echo "Starting processing for $files"
echo $outdir/$i.merged.vcf

# Index the BCF files
for j in $files
do nohup bcftools index $j &
done

numProc=$(ps auxw|grep -i bcftools|grep -v grep|wc -l)
while [ $numProc -gt 20 ]
do sleep 1
numProc=$(ps auxw|grep -i bcftools|grep -v grep|wc -l)
done


done
