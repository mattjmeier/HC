#!/bin/bash
# Takes a directory as input, and merges all the BCF files from separate samples, keeping each chromosome as a different file.
# This script assumes that the BCFs are indexed. See part 1 of the script to index.

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

numProc=$(ps auxw|grep -i bcftools|grep -v grep|wc -l)
while [ $numProc -gt 10 ]
do sleep 1
numProc=$(ps auxw|grep -i bcftools|grep -v grep|wc -l)
done

# Start a screen for each 
# screen -d -m -S "$i.$samples".merge_bcfs bash -c "bcftools merge $files -o $outdir/$i.merged.vcf"

##### PLEASE NOTE, I HAVE COMMENTED THIS OUT BECAUSE YOU MUST CAREFULLY SPECIFICY ARGUMENTS FOR MERGING ####
##### YOU SHOULD DECIDE IF EACH VARIANT GOES ON A SEPARATE LINE, OR IF SCRICT OUTPUT IS REQUIRED ###########
##### THERE ARE MANY OPTIONS SO PLEASE CHOOSE CAREFULLY ####################################################
# nohup bcftools merge $files -o $outdir/$i.merged.vcf &

done
