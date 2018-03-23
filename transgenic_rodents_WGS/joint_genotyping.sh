#!/bin/bash

# Runs joint genotyping for each chromosome


# Specify parent directory that joint genotyping should be done on
echo "On what directory should joint genotyping be run?"
read -e rootdir
echo "Where should the output go?"
read -e outdir

#for i in $(cat  /s0/ngs/references/mouse/ucsc/mm10/index/samtools/chr_list.txt)
#for i in $(cat  /s0/ngs/references/mouse/ucsc/mm10/index/samtools/chr_list.txt | grep -v "chr4$" | grep -v "chr7$" | grep -v "chr8$")
for i in $(echo chr13 )
#for i in $(echo chrM)

do files=$(find $rootdir -name "*$i.*.gz" | grep -vi old)

echo "Files found for each chromosome:"
echo $files 

infileArray=()
for j in $files
do
   echo $j
   filename=$(printf "%svariant %s " "--" "$j")
   infileArray+=($filename)
done
# printf '%s\n' "${infileArray[@]}"
variantFiles=$(echo "${infileArray[@]}")
echo $variantFiles

java -Xmx20G -jar $GATK36 \
   -T GenotypeGVCFs \
   -R $mm10 \
   -nt 32 \
   $variantFiles \
   -o $outdir/genotyped.$i.vcf.gz

done
