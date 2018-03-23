#!/bin/bash
# Merges multiple gVCFs using bcftools (experimental version)

echo "What is the name of the mother?"
read mother

echo "What is the name of the father?"
read father

echo "What is the name of the child?"
read child

mom_gvcf=$(find /s1/ngs/DNA/pedigrees/analysis/pedigree_BaP/tvc_output/Oct23/ -regex "\/s1\/ngs\/DNA\/pedigrees\/analysis\/pedigree_BaP\/tvc_output\/Oct23\/$mother\/TSVC_variants.genome.vcf.gz")
dad_gvcf=$(find /s1/ngs/DNA/pedigrees/analysis/pedigree_BaP/tvc_output/Oct23/ -regex "\/s1\/ngs\/DNA\/pedigrees\/analysis\/pedigree_BaP\/tvc_output\/Oct23\/$father\/TSVC_variants.genome.vcf.gz")
kid_gvcf=$(find . -regex ".*\.\/$child\/TSVC_variants.genome.vcf.gz")

echo "Mom:"
echo $mom_gvcf
echo "Dad:"
echo $dad_gvcf
echo "Kid:"
echo $kid_gvcf

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

echo "Merging gVCFs (bcftools)..."
$bcftoolsexp merge -O z -g $kid_gvcf $mom_gvcf $dad_gvcf > $child/$child.trio.merged.vcf.gz


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

