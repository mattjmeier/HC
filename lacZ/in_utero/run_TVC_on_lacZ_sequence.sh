#!/bin/bash

# RUN TORRENT VARIANT CALLER ON MULTIPLE BAM FILES, USING A SPECIFIC ALLELE FREQUENCY FOR EACH FILE

# Read in where the number of plaques for each sample can be found
#echo "Where are the plaque counts for this group stored?"
#read counts_file
counts_file=$1

# What directory filled with BAM files to run on
#echo "Where are the bam files you would like to run?"
#read bam_dir
bam_dir=$2

# Store BAM files in variable
files=$(find $bam_dir -name "IonXpress*.bam" | sort)
# Make an array for BAM files
files_array=($files)
for f in $files_array; do echo $f; done

# How many files are there?
echo "There are this many BAM files:"
numfiles=${#files_array[@]}
echo $numfiles


# Grab plaque count for file (all samples)
plaques=$(cat $counts_file | grep -Po [0-9]+$)

# Genotype frequency based on plaque count for each barcode; divided by two to reduce stringency.
gt=$(for i in $plaques; do echo "scale=3; 1/$i/2" | bc -l; done)
# Make an addressable array from the string
gt_array=($gt)

# Calculate appropriate coverage for each sample - 200X coverage per plaque
cov=$(for i in $plaques; do echo "scale=6; 200*$i" | bc -l; done)
# Make an addressable array from the string
cov_array=($cov)

# Grab sample name from counts file
sample=$(cat $counts_file | awk '{print $2}')
# Make an addressable array from the string
sample_array=($sample)


# Grab length of array
length=${#gt_array[@]}
echo "There are this many genotypes in the list"
echo $length

# Sanity check of array
for i in ${gt_array[@]}; do echo $i; done


for (( i=0; i<${length}; i++ ));

do

echo "Running for file number $i, using genotype frequency ${gt_array[$i]} on bamfile ${files_array[$i]}, outputting to file ${sample_array[$i]}"

# Run torrent variant caller using gt for each file
screen -d -m -S "$i.${sample_array[$i]}" bash -c "tvc -n 2 -r /home/ngs/references/lacz/lacZFlncCorr/ReferenceGenomeLacZ.fasta -t /home/ngs/references/lacz/lacZFlncCorr/laczGeneRegion.bed --downsample-to-coverage ${cov_array[$i]} --gen-min-alt-allele-freq ${gt_array[$i]} --gen-min-indel-alt-allele-freq ${gt_array[$i]} --snp-min-allele-freq ${gt_array[$i]} --indel-min-allele-freq ${gt_array[$i]} -b ${files_array[$i]} -o ${sample_array[$i]}.vcf 2>&1 | tee ${sample_array[$i]}.log; exit"

numProc=$(ps auxw|grep -i tvc|grep -v grep|wc -l)

while [ $numProc -gt 90 ]
do sleep 5
numProc=$(ps auxw|grep -i tvc|grep -v grep|wc -l)
done

done

# Pipe stdout to logfile named $fname.log
