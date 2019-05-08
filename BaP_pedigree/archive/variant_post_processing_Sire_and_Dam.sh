#!/bin/bash

echo "What is the sire of the family?"
read sire
echo "What is the dam of the family?"
read dam
echo "Where is your VCF file of mendelian violations?"
read -e vcf

# Make file for sire's BAM files
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$sire"_*.bam > /tmp/$sire.list
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$dam"_*.bam > /tmp/$dam.list

# Generate BED file:
cat $vcf | grep -v "^#" | awk ' { print $1"\t"$2-1"\t"$2 } ' > /tmp/$vcf.bed

# Generate alternate allele list:
# cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V > /tmp/$vcf.alts

# Get base info
java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o /tmp/$vcf.sire.allele_depth    -I /tmp/$sire.list --printBaseCounts    -L /tmp/$vcf.sire.bed
java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o /tmp/$vcf.dam.allele_depth    -I /tmp/$dam.list --printBaseCounts    -L /tmp/$vcf.dam.bed

# while read line; do cat $line | awk '{ print $1} '; done < ^CF1.trio.violations.filtered.vcf > /tmp/file

# For each line in the BED region file
# Search for coordinates in 

cat /tmp/$vcf.sire.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $vcf.sire.basecounts
cat /tmp/$vcf.dam.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $vcf.dam.basecounts

# Parse the allele depth (basecounts) output from GATK; look at 
cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $vcf.sire.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$vcf.sire.annotations
cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $vcf.dam.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$vcf.dam.annotations



# Add to above code to print out specific base counts, for troubleshooting.
# base "\tA"$3"\tC"$4"\tG"$5"\tT"$6


# FORMAT FOR ANNOTATIONS FILE
# CHR     FROM   TO      ANNOTATION 
# 1        12345  22345   gene1 
# 1        67890  77890   gene2 

bgzip /tmp/$vcf.sire.annotations
tabix -s 1 -b 2 -e 3 /tmp/$vcf.sire.annotations.gz

bgzip /tmp/$vcf.dam.annotations
tabix -s 1 -b 2 -e 3 /tmp/$vcf.dam.annotations.gz


cat $vcf | vcf-annotate -a /tmp/$vcf.sire.annotations.gz \
   -d key=INFO,ID=SIRE_AF_FILTER,Number=1,Type=Integer,Description='This means the Sire has been found to have 0 contaminating bases of the de novo allele, when set to 0.' \
   -c CHROM,FROM,TO,INFO/SIRE_AF_FILTER > $vcf.Sire_AF_flagged.vcf


cat $vcf | vcf-annotate -a /tmp/$vcf.dam.annotations.gz \
   -d key=INFO,ID=DAM_AF_FILTER,Number=1,Type=Integer,Description='This means the Dam has been found to have 0 contaminating bases of the de novo allele, when set to 0.' \
   -c CHROM,FROM,TO,INFO/DAM_AF_FILTER > $vcf.Sire_and_Dam_AF_flagged.vcf


# if ( type is not del && alt basecounts are > 0 ) then annotate the VCF Filter column with FAIL


#!/bin/bash

# Provide a VCF file
# Provide list of associated BAM files
# Gets positions
# Returns bases

#echo "Where is your list of BAM files?"
#read -e bamfiles
#echo "Where is the VCF to process?"
#read -e vcf

# Get positions of de novo SNPs
#cat $vcf | grep -v '^#' | awk '{print $1":"$2"-"$2}' > $vcf.intervals
# Get positions and alt allele
#cat 32CF1.vcf | grep -v '^#' | awk '{print $1":"$2"\t"$5}' | less
# GRAB FILTERED POSITIONS
#vcftools --vcf 32CF1.vcf --bed passed_base_cutoffs.bed --out 32CF1_base_cutoff --recode --keep-INFO-all

# Get pileup stats for each position
#while read pos; do
#~/DNA/pedigree_BaP/bcf_by_sample/dng_ouput/vcfs-merged/piledriver/piledriver/bin/bamtools piledriver -list $bamfiles -region $pos
#done < $vcf.regions > $vcf.pileup

