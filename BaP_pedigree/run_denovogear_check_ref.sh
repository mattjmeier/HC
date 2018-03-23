#!/bin/bash

# get working directory
echo "What directory are the bcfs in?"
read -e rootdir

# get files that belong to same family
echo "What sample names should be processed? Please provide a list of exact names, separated by the pipe symbol, |"
read samples
echo "What should we name this family?"
read family

# get pedigree file
echo "What pedigree file should be used for defining the trio relationships?"
read -e pedigree

# get output directory
echo "Where should the output go?"
read -e outdir


# loop through chromosomes
chromosomes=$(find $rootdir -name "*.bcf" | sed 's/^.*\(chr[0-9A-Za-z]*\)\.bcf/\1/g' | sort | uniq)

for i in $chromosomes
do
echo $i

# files=$(find . -name "*.$i.bcf")/DNA/pedigree_BaP/bcf_by_sample
files=$(find $rootdir -maxdepth 1 -regextype posix-egrep -regex ".*($samples).$i.bcf") # This one could be used to find for specific samples...
myArray=()
for j in $files
do
myArray+=($j)
done

echo "Starting processing for ${myArray[@]}"
echo $outdir $family $i $pedigree

#screen -d -m -S $family.$i.denovogear bash -c "Sire32|32CF1|Dam32C
#sleep 10

nohup $bcftoolsexp merge -m none ${myArray[@]} -Ou | $bcftoolsexp norm -m +any --check-ref w -f $mm10min | /utils/appz/denovogear/denovogear-v1.1.1-Linux-x86_64/bin/dng dnm auto --ped $pedigree --output_vcf  $outdir/$family.$i.dng-output.vcf --vcf - > $outdir/$family.$i.dng-output.txt &

#"

# Dont start too many DNG instances
numProc=$(ps auxw|grep -i dng|grep -v grep|wc -l)
while [ $numProc -gt 20 ]
do sleep 1
numProc=$(ps auxw|grep -i dng|grep -v grep|wc -l)
done


done

