#!/bin/bash
# Takes a directory as input, and merges all the BCF files from separate samples, keeping each chromosome as a different file.
# Outputs uncompressed VCF files for input to DeNovoGear

echo "What directory (or directories, separated by spaces) are your VCF files, split by chromosome, located in?"
read -e rootdir

echo "What sample names should be processed? Please provide a list of exact names, separated by the pipe symbol, |"
read samples

echo "Where should the merged VCF files be written?"
read -e outdir

chromosomes=$(find $rootdir -maxdepth 1 -name "*.vcf" | sed 's/.*\(chr[0-9A-Za-z]*\).*.vcf/\1/g' | sort | uniq)


#### LOOP THROUGH EACH CHROMOSOME ####

for i in $chromosomes
do
echo "Chromosome: $i"

# Find files for all samples to be merged, then print them in proper format for GATK
infileArray=()
files=$(find $rootdir -maxdepth 1 -regextype posix-egrep -regex ".*($samples).*$i.*.vcf")
for j in $files; do echo $j; infileArray+=($j); done
vcfs=$(printf " -V %s" "${infileArray[@]}")

#echo "Starting processing for $vcfs"
echo $outdir/$i.merged.vcf

numProc=$(ps auxw|grep -i CombineVariants|grep -v grep|wc -l)
while [ $numProc -gt 15 ]
do sleep 1
numProc=$(ps auxw|grep -i CombineVariants|grep -v grep|wc -l)
done

echo $i
echo $samples
echo $vcfs
# Start a screen for all 
screen -d -m -S "$i"."$samples".mergeVCFs bash -c "

java -jar $GATK \
   -T CombineVariants \
   -R $mm10min \
   $vcfs \
   -o $outdir/$i.merged.vcf \
   -genotypeMergeOptions UNSORTED

"


done
