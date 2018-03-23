#!/bin/bash

# Takes a directory full of recalibrated BAM files
# Runs GATK Haplotype Caller in DISCOVERY MODE
# COMPLETE FOR EVERY SAMPLE

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."

echo "What directory would you like to run GATK Variant Discovery on? This should be done on a per-sample basis."
#read -e rootdir
rootdir=/s1/ngs/DNA/pedigrees/
echo "What sample would you like to process?"
read -e sample
echo "What directory would you like to put the output files in? Please provide the full path."
#read -e outdir
outdir=/s1/ngs/DNA/pedigrees/gatk_variant_discovery

files=$(find $rootdir -maxdepth 1 -name "*.recal_reads.bam" -type f | grep -v old | sort)


# Count number of files in input; report them as raw filenames
num=$(for i in $files; do echo $i; done | wc -l)
echo "Selected input directory has the following $num files:"
# Create an array of input files
infileArray=()
for i in $files
do
echo $i
# Add inputFile to array
sm=$(samtools view -h $i | head -n 100 | sed -n -e 's/^.*SM:\([0-9A-Za-z:]*\).*$/\1/p' | uniq)
# If BAM file sample name matches the requested sample name, add to array
if [ "$sm" == "$sample" ]; then
infileArray+=($i)
fi
done
# Report the input file array
echo "The input directory has been parsed to the following ${#infileArray[@]} files:"

printf '%s\n' "${infileArray[@]}"

# Store header from BAM file
# NOTE they must be aligned to the same reference sequence, because we are only taking the header from ONE
echo "Grabbing header from ${infileArray[0]}"
bamHeader=$(samtools view -H ${infileArray[0]})

# Obtain a list of all the sequences from the BAM headers
# Grab lines starting with @SQ, then take the second column and used sed to remove the SN: part
myList=$(echo "$bamHeader" | awk -F" " '$1 == "@SQ"{print $2}' | sed 's/SN\://g')

echo "The IDs are:"
for id in $myList
do
echo $id
done

read -r -p "Please ensure the above is correct before you continue!!! [y/N] " response
case $response in
    [yY][eE][sS]|[yY]) 

        echo "Going ahead with analysis..."
	# You have 3 seconds to kill it
	sleep 3

###########################################################
###########################################################
###  THIS IS THE CODE TO RUN GATK #########################
###########################################################
###########################################################

for i in $myList

do
echo "Starting job for $sample, $i..."

# Store file names to be processed from infileArray as a list of -I inputs for GATK
files=$(printf " -I %s" "${infileArray[@]}")

# Start a screen for each chromosome
screen -d -m -S "$sample.$i".Haplotype.Caller.discovery bash -c "

java -Xmx20g -jar $GATK3446 \
    -T HaplotypeCaller \
    -R $mm10min \
    $files \
    --emitRefConfidence GVCF \
    -L $i \
    -o $outdir/$sample.$i.g.vcf
"

# A few extra options for haplotype caller that might be useful
#    -stand_emit_conf 10 \
#    -stand_call_conf 30 \
#    --sample_name $sample \
#    --dbsnp $mm10_snps \
#    --genotyping_mode DISCOVERY \

# Make sure things start ok
# sleep 10

# Don't start more than 10 instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 10 ]
do sleep 1
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
done

done

###########################################################
###########################################################
###  END GATK #############################################
###########################################################
###########################################################

        ;;
    *)
        echo "Quitting now"
        ;;
esac
