#!/bin/bash

# Merges multiple gVCFs using bcftools (experimental version)

echo "What is the name of the mother?"
read mother

echo "What is the name of the father?"
read father

echo "What is the name of the child?"
read child

mom_gvcf=$(find . -regex ".*\.\/$mother\/TSVC_variants.genome.vcf.gz")
dad_gvcf=$(find . -regex ".*\.\/$father\/TSVC_variants.genome.vcf.gz")
kid_gvcf=$(find . -regex ".*\.\/$child\/TSVC_variants.genome.vcf.gz")

echo "Mom:"
echo $mom_gvcf
echo "Dad:"
echo $dad_gvcf
echo "Kid:"
echo $kid_gvcf

sire=$father
dam=$mother

# Make file for sire's BAM files
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$sire"_*.bam > /tmp/$sire.list
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$dam"_*.bam > /tmp/$dam.list

echo "The list of BAM files associated with the sire is:"
cat /tmp/$sire.list
echo "The list of BAM files associated with the dam is:"
cat /tmp/$dam.list

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


$bcftoolsexp merge -O z -g $kid_gvcf $mom_gvcf $dad_gvcf | tee $child/$child.trio.merged.vcf.gz | $bcftoolsexp view | java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )" > $child/$child.trio.violations.vcf


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

cat /tmp/$vcf.sire.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$vcf.sire.basecounts
cat /tmp/$vcf.dam.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$vcf.dam.basecounts

# Parse the allele depth (basecounts) output from GATK; look at 
cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$vcf.sire.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$vcf.sire.annotations
cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$vcf.dam.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$vcf.dam.annotations



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
   -c CHROM,FROM,TO,INFO/SIRE_AF_FILTER > $child/$vcf.Sire_AF_flagged.vcf


cat $vcf | vcf-annotate -a /tmp/$vcf.dam.annotations.gz \
   -d key=INFO,ID=DAM_AF_FILTER,Number=1,Type=Integer,Description='This means the Dam has been found to have 0 contaminating bases of the de novo allele, when set to 0.' \
   -c CHROM,FROM,TO,INFO/DAM_AF_FILTER > $child/$vcf.Dam_AF_flagged.vcf



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

