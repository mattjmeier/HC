#!/bin/bash

# Takes MVs and filters them

rootdir=~/DNA/pedigree_BaP/tvc_output/2017_additional_offspring

echo "What is the name of the dam?"
read dam

echo "What is the name of the sire?"
read sire

echo "What is the name of the offspring?"
read child

# Make file for BAM lists
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$sire"_*.bam > $rootdir/$child/$sire.list
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$dam"_*.bam > $rootdir/$child/$dam.list
find /s2/ngs/DNA/pedigrees/BaP/raw/ -name *_"$child"_*.bam >  $rootdir/$child/$child.list


echo "The following BAM files have been identified for the sire:"
cat $rootdir/$child/$sire.list
cat $rootdir/$child/$sire.list   > $rootdir/$child/$child.allbams.list
echo "The following BAM files have been identified for the dam:"
cat $rootdir/$child/$dam.list 
cat $rootdir/$child/$dam.list   >> $rootdir/$child/$child.allbams.list
echo "The following BAM files have been identified for the child:"
cat $rootdir/$child/$child.list
cat $rootdir/$child/$child.list >> $rootdir/$child/$child.allbams.list
echo "The combined files are:"
cat $rootdir/$child/$child.allbams.list

read -r -p "Please ensure the above is correct before you continue!!! [y/N] " response
case $response in
    [yY][eE][sS]|[yY])

        echo "Going ahead with analysis..."
        # You have 1 second to kill it
        sleep 1

###########################################################
###########################################################
###  START COMMANDS #######################################
###########################################################
###########################################################




# Filter based on various metrics - depth, quality, strand bias, etc.
# Remove variants in segmentally duplicated regions

echo "Removing variants found in segmentally duplicated regions and filtering based on hard-coded cutoff values for depth, quality, strand bias, and allele fraction..."
cat $rootdir/$child/$child.violations.raw.vcf | \
java -jar /utils/appz/snpEff/SnpSift.jar filter "  \
	isVariant( GEN[$child] ) && isRef( GEN[$dam] ) && isRef( GEN[$sire] ) && isHom( GEN[$dam] ) && isHom( GEN[$sire] ) && \
	(QUAL >= 10) &&  (MIN_DP >= 10) &&  (DP >=10) &&  (STBP >= 0.05) && (QD >= 5) && (GEN[$child].AF > 0.25 ) " | \
java -jar /utils/appz/snpEff/SnpSift.jar intervals -x /s0/ngs/references/mouse/ucsc/mm10/index/samtools/SegmentalDups_mm10.bed > $rootdir/$child/$child.violations.filtered.vcf

echo "Done."

# This should pare down the initial raw MV file to one containing about 5% of the original number.  (Approx 400,000 down to 10,000)

# Next we will filter based on parental alleles, to remove sites where the variant is observed in the parents.

# Generate a BED file for all potential de novo sites for the family:

echo "Making a BED file with variant positions from filtered Mendelian violations file..."
echo "Insertions..."
vcf2bed -t  < $rootdir/$child/$child.violations.filtered.vcf >		$rootdir/$child/$child.filtered_violations_ins.bed
echo "Deletions..."
vcf2bed -n  < $rootdir/$child/$child.violations.filtered.vcf >		$rootdir/$child/$child.filtered_violations_del.bed
echo "SNVs..."
vcf2bed -v  < $rootdir/$child/$child.violations.filtered.vcf >		$rootdir/$child/$child.filtered_violations_snv.bed
echo "Combining..."
cat $rootdir/$child/$child.filtered_violations_*.bed 			| sortBed >					$rootdir/$child/$child.filtered_violations.combined.sorted.bed
cat $rootdir/$child/$child.filtered_violations.combined.sorted.bed 	| awk ' {print $1"\t"$2"\t"$3} ' >		$rootdir/$child/$child.filtered_violations.trunc.bed
echo "Done."

# Alternate method for making bed file (do not use for deletions)
# cat $rootdir/$child/$child.violations.filtered.vcf | grep -v "^#" | awk ' { print $1"\t"$2-1"\t"$2 } ' > $rootdir/$child/$child.bed

# Get base info from BAM files using GATK
echo "Running GATK DepthOfCoverage..."
while true
do
nohup  java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o $rootdir/$child/$child.sire.raw.basecounts     -I $rootdir/$child/$sire.list  --printBaseCounts   -L $rootdir/$child/$child.filtered_violations.trunc.bed &
        pid1=$!
nohup  java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o $rootdir/$child/$child.dam.raw.basecounts      -I $rootdir/$child/$dam.list   --printBaseCounts   -L $rootdir/$child/$child.filtered_violations.trunc.bed &
        pid2=$!
nohup  java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o $rootdir/$child/$child.child.raw.basecounts    -I $rootdir/$child/$child.list --printBaseCounts   -L $rootdir/$child/$child.filtered_violations.trunc.bed &
        pid3=$!
        wait $pid1
        wait $pid2
        wait $pid3
done
echo "Done."

# Get base info from BAM files using samtools
echo "Running samtools mpileup..."
	samtools mpileup -f $mm10min  -b $rootdir/$child/$child.allbams.list -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -uvl $rootdir/$child/$child.filtered_violations.trunc.bed > $rootdir/$child/$child.samtools.vcf
echo "Done."


# while read line; do cat $line | awk '{ print $1} '; done < 32CF1.trio.violations.filtered.vcf > $rootdir/$child/file

# For each line in the BED region file
echo "Making basecounts files..."

while true
do
        cat $rootdir/$child/$child.sire.raw.basecounts | grep "^chr" | \
        	awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.sire.basecounts
        pid1=$!
        cat $rootdir/$child/$child.dam.raw.basecounts | grep "^chr" | \
        	awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.dam.basecounts
        pid2=$!
        cat $rootdir/$child/$child.child.raw.basecounts | grep "^chr" | \
        	awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.child.basecounts
        pid3=$!
        wait $pid1
        wait $pid2
        wait $pid3
done
echo "Done."

###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################

# Parse the allele depth (basecounts) output from GATK

#echo "Making annotations for the sire allele frequency..."
#cat $rootdir/$child/$child.violations.filtered.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.sire.basecounts - | 		  \
#	awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; \
#	if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > $rootdir/$child/$child.sire.annotations
	
#echo "Making annotations for the dam allele frequency..."
#cat $rootdir/$child/$child.violations.filtered.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.dam.basecounts - | 		  \
#	awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; \
#	if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > $rootdir/$child/$child.dam.annotations
	
#echo "Making annotations for the child allele frequency..."
#cat $rootdir/$child/$child.violations.filtered.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.child.basecounts - | 	  \
#	awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; \
#	if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > $rootdir/$child/$child.child.annotations

# Add to above code to print out specific base counts, for troubleshooting.
# base "\tA"$3"\tC"$4"\tG"$5"\tT"$6

###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################

# FORMAT FOR ANNOTATIONS FILE
# CHR     FROM   TO      ANNOTATION
# 1        12345  22345   gene1
# 1        67890  77890   gene2

###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################


echo "Running bgzip for sire annotations file..."
bgzip $rootdir/$child/$child.sire.annotations
echo "Running tabix for sire annotations file..."
tabix -fs 1 -b 2 -e 3 $rootdir/$child/$child.sire.annotations.gz

echo "Running bgzip for dam annotations file..."
bgzip $rootdir/$child/$child.dam.annotations
echo "Running tabix for dam annotations file..."
tabix -fs 1 -b 2 -e 3 $rootdir/$child/$child.dam.annotations.gz

echo "Applying annotations file to $child_gvcf..."
cat $rootdir/$child/$child.violations.filtered.vcf | vcf-annotate -a $rootdir/$child/$child.sire.annotations.gz 				\
   -d key=INFO,ID=SIRE_ALT_READS,Number=1,Type=Integer,Description='Number of reads from the Sire containing the de novo alternate allele' 	\
   -c CHROM,FROM,TO,INFO/SIRE_ALT_READS | vcf-annotate -a $rootdir/$child/$child.dam.annotations.gz 						\
   -d key=INFO,ID=DAM_ALT_READS,Number=1,Type=Integer,Description='Number of reads from the Dam containing the de novo alternate allele'   	\
   -c CHROM,FROM,TO,INFO/DAM_ALT_READS > $child/$child.Sire_and_Dam_AF_flagged.vcf

#   > $child/$child.Sire_AF_flagged.vcf

# cat $child/$child.Sire_AF_flagged.vcf | vcf-annotate -a $rootdir/$child/$child.dam.annotations.gz \
#   -d key=INFO,ID=DAM_ALT_READS,Number=1,Type=Integer,Description='This means the Dam has been found to have 0 contaminating bases of the de novo allele, when set to 0.' \
#   -c CHROM,FROM,TO,INFO/DAM_ALT_READS > $child/$child.Sire_and_Dam_AF_flagged.vcf

echo "Making a summary file of filtered Mendelian violations..."
cat $child/$child.Sire_and_Dam_AF_flagged.vcf | grep -v "\./\." | java -jar /utils/appz/snpEff/SnpSift.jar filter " (DAM_ALT_READS == 0) & (SIRE_ALT_READS == 0) "  | 				\
	tee $child/$child.MVs.NO_MISSING.vcf | java -jar /utils/appz/snpEff/SnpSift.jar filter " ( QUAL > 10 ) & ( MIN_DP > 10 ) " | tee $child/$child.MVs.QUAL_10.MIN_DP_10.NO_MISSING.vcf | 	\
	grep -v "^#" | awk '{print $1","$2}' > $child/$child.MVs.QUAL_10.MIN_DP_10.NO_MISSING.csv


###########################################################
###########################################################
###  END COMMANDS #########################################
###########################################################
###########################################################

        ;;
    *)
        echo "Quitting now"
        ;;
esac

