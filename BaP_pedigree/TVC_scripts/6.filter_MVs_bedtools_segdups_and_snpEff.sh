#!/bin/bash

# Takes MVs and filters them

rootdir=~/DNA/pedigree_BaP/tvc_output/2017_additional_offspring

while read dam sire child; do

echo "Mom: $dam"
echo "Dad: $sire"
echo "Kid: $child"

echo "File to be processed:"
find $rootdir -name "*$child.violations.raw.vcf"
echo "Should match:"
echo $rootdir/$child/$child.violations.raw.vcf

echo "Sleeping..."
sleep 5
echo "Done sleeping."

###########################################################
###########################################################
###  START COMMANDS #######################################
###########################################################
###########################################################

# Filter based on various metrics - depth, quality, strand bias, etc.
# Remove variants in segmentally duplicated regions

echo "Removing variants found in segmentally duplicated regions..."

cat 		     $rootdir/$child/$child.violations.raw.vcf 	 | grep '^#' 									>  $rootdir/$child/$child.violations.segdups_filtered.vcf
bedtools subtract -a $rootdir/$child/$child.violations.raw.vcf   -b /s0/ngs/references/mouse/ucsc/mm10/index/samtools/SegmentalDups_mm10.bed  	>> $rootdir/$child/$child.violations.segdups_filtered.vcf

echo "Done."



echo "Filtering based on hard-coded cutoff values for depth, quality, and offspring allele fraction..."
# Add strand bias filter for STBP?

#cat $rootdir/$child/$child.violations.segdups_filtered.vcf | java -jar /utils/appz/snpEff/SnpSift.jar filter  \
#	" 	isVariant( GEN[$child] ) && isRef( GEN[$dam] ) && isRef( GEN[$sire] ) && isHom( GEN[$dam] ) && isHom( GEN[$sire] ) && \
#	(QUAL >= 10) &&  (MIN_DP >= 10) &&  (DP >=10) &&  (QD >= 5) && (GEN[$child].AF > 0.25 ) 	"			>  $rootdir/$child/$child.violations.segdups_and_hard_filters.vcf

vcf_file="${rootdir}/${child}/${child}.violations.segdups_filtered.vcf"


echo "VCF file input: echo ${vcf_file}"
echo "output: $rootdir/$child/$child.violations.segdups_and_hard_filters.vcf"
echo "offspring: $child"
echo "dam: $dam"
echo "sire: $sire"

cat $vcf_file | java -jar /utils/appz/snpEff/SnpSift.jar filter  " (QUAL >= 10) && (MIN_DP >= 10) && (DP >=10) && (QD >= 5) && isVariant( GEN[$child] ) && isRef( GEN[$dam] ) && isRef( GEN[$sire] ) && isHom( GEN[$dam] ) && isHom( GEN[$sire] ) "	>  $rootdir/$child/$child.violations.segdups_and_hard_filters.vcf


echo "Done."

# This should pare down the initial raw MV file to one containing about 5% of the original number.  (Approx 400,000 down to 10,000)

###########################################################
###########################################################
###  END COMMANDS #########################################
###########################################################
###########################################################

done < families.txt


