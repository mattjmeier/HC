#!/bin/bash

# Takes a directory of VCF files (each containing at least 3 biological samples)
# Runs BCF tools to call SNPS
# Runs GATK PhaseByTransmission to determine family characteristics

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "This script will take a directory of VCF files containing RAW calls, preferably split by chromosome with each containing SNP information for at least one trio.  It will run GATK PhaseByTransmission, optionally outputting a mendelian violations file to look for de novo mutations."
printf '\n\n'
echo "REMOVE -maxdepth 1 ARGUMENT FROM FIND TO DESCEND INTO DIRECTORIES, IF DESIRED.  Currently, the check for processed files parses filenames to the first period; ensure there are no abmiguous names in your directory before proceeding."

echo "What directory would you like to process? IT SHOULD CONTAIN VCFS OF RAW READS."
read -e rootdir
echo "What directory would you like to put the output files in?"
read -e outdir
echo "What pedigree file should be used for defining the trio relationships?"
read -e pedigree
echo "Would you like to grep for a particular pattern? If so, enter it; if not, leave blank and press enter"
read pattern

files=$(find $rootdir -maxdepth 1 -regex ".*.[bv]cf" -type f | grep $pattern) # | grep -v old | grep -v analysis_files | grep $pattern | sort)
processed=$(find $outdir -maxdepth 1 -name "*.filtered.vcf" -type f | sort)

num=$(for i in $files; do echo $i; done | wc -l)
echo "Selected input directory has the following $num files:"
# Create an array of input files
infileArray=()
for i in $files
do
echo $i
# Substitute for base filename
inputFile=${i##*/}
inputFile=${inputFile%%.*}
# Add inputFile to array
infileArray+=($inputFile)
done
# Report the input file array
echo "The input directory has been parsed to the following names:"
echo ${infileArray[@]}

# Count number of files that have been processed; report them as raw filenames
num=$(for i in $processed; do echo $i; done | wc -l)
echo "Files already processed in your selected output directory, $num files:"
doneArray=()
for i in $processed
do
echo $i
# Substitute for base filename
doneFile=${i##*/}
doneFile=${doneFile%%.*}
# Add doneFile to array
doneArray+=($doneFile)
done
# Report the processed file array
echo "List of processed files:"
echo ${doneArray[@]}

echo "Output from printf of infileArray:"
printf '%s\n' "${infileArray[@]}"
echo "Output from printf of doneArray:"
printf '%s\n' "${doneArray[@]}"

# Takes the files unique to the list of input files
# Ignores files that are present in both lists
# Prepends the input directory and appends the extension to make sure the runlist points to the correct files
runlist=()
runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) )

FILEPATHS=()
for i in $runlist
do
j=$(find $rootdir -maxdepth 1 -regex ".*$i.*.[bv]cf" -type f | grep $pattern ) # | grep -v old | grep -v analysis_files)
FILEPATHS+=($j)
done

num=$(for i in "${FILEPATHS[@]}"; do echo $i; done | wc -l)
echo "The final FILEPATHS list has the following $num files:"
for i in "${FILEPATHS[@]}"
do
echo $i
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

for i in "${FILEPATHS[@]}"

do

fname=${i##*/}
fname=${fname%%.vcf*}
echo "Now starting screen for sample $fname, file path is $i"
echo "The input directory is $rootdir"
echo "The output directory is $outdir"

# Start a screen session detatched and named by the file being processed
# Run commands in new bash shell

# OPTIONS PRESET FOR ANALYSIS USING TRIOCALLER
#screen -d -m -S "$fname".BCFtools bash -c "
#bcftools view --exclude-types indels,mnps,other $i | bcftools call -v --constrain trio -S $pedigree -m - | bcftools view --exclude-uncalled - | vcftools --max-missing 1.0 --vcf - --recode --out $outdir/$fname.raw_calls.vcf
#"


# OPTIONS PRESET FOR ANALYSIS USING GATK PhaseByTransmission PIPELINE
# Currently pointing to the devel version of BCFTools because of a bug in the working version that results in badly formatted VCF files - sometimes outputs ./. when it should be 0/0 in GT field.

#screen -d -m -S "$fname".BCFtools.PBT bash -c "

#sleep 10

java -jar $GATK3446 \
 -T PhaseByTransmission \
 -R $mm10min \
 -V $i \
 -ped $pedigree \
 -mvf $outdir/$fname.violations.txt \
 -l DEBUG \
 -o $outdir/$fname.phased.vcf

java -jar $GATK \
   -T SelectVariants \
   -R $mm10min  \
   -V $outdir/$fname.phased.vcf \
   -ped $pedigree \
   -mv -mvq 20 \
   -o $outdir/$fname.violations_filtered.vcf

#"

# USE THIS TO ADD THE MISSING GQ FIELD FROM THE HEADER... I HAVE NO IDEA WHY THIS ISN'T ADDED PROPERLY.
# for i in $files; do sed -i.bak '/Description="Genotype"/a ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality"> ' $i; done

# ONE LINER TO RUN PHASE BY TRANSMISSION
# for i in $files; do java -jar $GATK  -T PhaseByTransmission  -R $mm10min  -V $i  -ped $pedigree -mvf $i.violations.txt -o $i.phased.vcf ; done

# Make sure things start ok
sleep 5

# Don't start more than X instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 10 ]
do sleep 5
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
