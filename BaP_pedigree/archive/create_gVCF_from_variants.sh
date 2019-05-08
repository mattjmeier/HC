#!/bin/bash

# Takes a list of sample names

# For each sample name, 1) grab BAM file list, and 2) merge into gVCF

# Root directory for project BAM files
rootdir=/s0/ngs/DNA/pedigrees/BaP/raw/

echo "Where is the list of samples you would like to make gVCFs for?"
read -e sample_list

topdir=`pwd`

for i in $(cat $sample_list)

do

cd $topdir

path=$(find `pwd` -type d -name $i)

files=$(find $rootdir -name *_"$i"_*.bam | tr '\n' ' ' | rev | cut -c 2- | rev)

cd $path

nohup samtools depth $files | tvcutils unify_vcf --novel-tvc-vcf small_variants.vcf --novel-assembly-vcf indel_assembly.vcf --output-vcf TSVC_variants.vcf.gz --reference-fasta $mm10min --tvc-metrics tvc_metrics.json --input-depth stdin --min-depth  5 &

done
