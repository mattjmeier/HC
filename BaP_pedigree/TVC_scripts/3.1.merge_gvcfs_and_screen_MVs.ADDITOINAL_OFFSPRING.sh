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
kid=$child

# Make file for BAM lists
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$sire"_*.bam > /tmp/$sire.list
find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$dam"_*.bam > /tmp/$dam.list
#find /s0/ngs/DNA/pedigrees/BaP/raw/ -name *_"$kid"_*.bam > /tmp/$kid.list
find /s2/ngs/DNA/pedigrees/BaP/raw/ -name *_"$child"_*.bam >  /tmp/$kid.list

echo "The list of BAM files associated with the sire is:"
cat /tmp/$sire.list
echo "The list of BAM files associated with the dam is:"
cat /tmp/$dam.list
echo "The list of BAM files associated with the kid is:"
cat /tmp/$kid.list



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

echo "Merging gVCFs (bcftools) and calling Mendelian violations where kid is a variant, and parents are both hom-ref..."
$bcftoolsexp merge -O z -g $kid_gvcf $mom_gvcf $dad_gvcf | tee $child/$child.trio.merged.vcf.gz | $bcftoolsexp view | java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )" > $child/$child.trio.violations.vcf
echo "Done."

# Generate BED file:
echo "Making a BED file with variant positions from initial Mendelian violations file..."
cat $child/$child.trio.violations.vcf | grep -v "^#" | awk ' { print $1"\t"$2-1"\t"$2 } ' > /tmp/$child.bed
echo "Done."

# Generate alternate allele list:
# cat $vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V > /tmp/$vcf.alts


# Get base info
echo "Running GATK DepthOfCoverage..."
while true
do
	java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o /tmp/$child.sire.allele_depth    -I /tmp/$sire.list --printBaseCounts    -L /tmp/$child.bed
	pid1=$!
	java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o /tmp/$child.dam.allele_depth    -I /tmp/$dam.list --printBaseCounts    -L /tmp/$child.bed
	pid2=$!
	java -jar $GATK3446    -T DepthOfCoverage    -R $mm10min    -o /tmp/$child.kid.allele_depth    -I /tmp/$kid.list --printBaseCounts    -L /tmp/$child.bed
	pid3=$!
	wait $pid1
	wait $pid2	
	wait $pid3	
done



# while read line; do cat $line | awk '{ print $1} '; done < ^CF1.trio.violations.filtered.vcf > /tmp/file

# For each line in the BED region file
echo "Making basecounts files..."

while true
do
	cat /tmp/$child.sire.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.sire.basecounts
	pid1=$!
	cat /tmp/$child.dam.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.dam.basecounts
	pid2=$!
	cat /tmp/$child.kid.allele_depth | grep "^chr" | awk '{ split($1,POS,":"); split($5,A,":"); split($6,C,":"); split($7,G,":"); split($8,T,":"); print POS[1]":"POS[2]"\t"$4"\t"A[2]"\t"C[2]"\t"G[2]"\t"T[2] }' | sort -V > $child/$child.kid.basecounts
	pid3=$!
	wait $pid1
	wait $pid2
	wait $pid3
done
# Parse the allele depth (basecounts) output from GATK
echo "Making annotations for the sire allele frequency..."
cat $child/$child.trio.violations.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.sire.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$child.sire.annotations
echo "Making annotations for the dam allele frequency..."
cat $child/$child.trio.violations.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.dam.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$child.dam.annotations
echo "Making annotations for the kid allele frequency..."
cat $child/$child.trio.violations.vcf | grep -v "^#" | awk ' { print $1":"$2"\t"$5 } ' | sort -V | join $child/$child.kid.basecounts - | awk ' { if ($7 == "A") base=$3; else if ($7 == "C") base=$4; else if ($7 == "G") base=$5; else if ($7 == "T") base=$6; else base=" "; split($1,POS,":") ; if (base == 0 || base == 1) print POS[1]"\t"POS[2]"\t"POS[2]"\t"base } ' > /tmp/$child.kid.annotations


# Add to above code to print out specific base counts, for troubleshooting.
# base "\tA"$3"\tC"$4"\tG"$5"\tT"$6


# FORMAT FOR ANNOTATIONS FILE
# CHR     FROM   TO      ANNOTATION 
# 1        12345  22345   gene1 
# 1        67890  77890   gene2 
echo "Running bgzip for sire annotations file..."
bgzip /tmp/$child.sire.annotations
echo "Running tabix for sire annotations file..."
tabix -fs 1 -b 2 -e 3 /tmp/$child.sire.annotations.gz

echo "Running bgzip for dam annotations file..."
bgzip /tmp/$child.dam.annotations
echo "Running tabix for dam annotations file..."
tabix -fs 1 -b 2 -e 3 /tmp/$child.dam.annotations.gz

echo "Applying sire annotations file to $kid_gvcf..."
cat $child/$child.trio.violations.vcf | vcf-annotate -a /tmp/$child.sire.annotations.gz \
   -d key=INFO,ID=SIRE_AF_FILTER,Number=1,Type=Integer,Description='This means the Sire has been found to have 0 contaminating bases of the de novo allele, when set to 0; 1 contaminating base when set to 1; etc.' \
   -c CHROM,FROM,TO,INFO/SIRE_AF_FILTER | vcf-annotate -a /tmp/$child.dam.annotations.gz \
   -d key=INFO,ID=DAM_AF_FILTER,Number=1,Type=Integer,Description='This means the Dam has been found to have 0 contaminating bases of the de novo allele, when set to 0; 1 contaminating base when set to 1; etc.' \
   -c CHROM,FROM,TO,INFO/DAM_AF_FILTER > $child/$child.Sire_and_Dam_AF_flagged.vcf
   
   
#   > $child/$child.Sire_AF_flagged.vcf

# cat $child/$child.Sire_AF_flagged.vcf | vcf-annotate -a /tmp/$child.dam.annotations.gz \
#   -d key=INFO,ID=DAM_AF_FILTER,Number=1,Type=Integer,Description='This means the Dam has been found to have 0 contaminating bases of the de novo allele, when set to 0.' \
#   -c CHROM,FROM,TO,INFO/DAM_AF_FILTER > $child/$child.Sire_and_Dam_AF_flagged.vcf

echo "Making a summary file of filtered Mendelian violations..."
cat $child/$child.Sire_and_Dam_AF_flagged.vcf | grep -v "\./\." | java -jar /utils/appz/snpEff/SnpSift.jar filter " (DAM_AF_FILTER == 0) & (SIRE_AF_FILTER == 0) "  \
 | tee $child/$child.MVs.NO_MISSING.vcf | java -jar /utils/appz/snpEff/SnpSift.jar filter " ( QUAL > 10 ) & ( MIN_DP > 10 ) " | tee $child/$child.MVs.QUAL_10.MIN_DP_10.NO_MISSING.vcf | grep -v "^#" | awk '{print $1","$2}' > $child/$child.MVs.QUAL_10.MIN_DP_10.NO_MISSING.csv


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

