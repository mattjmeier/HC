#!/bin/bash

# Takes a directory full of realigned BAM files
# Runs GATK Base Recalibration
# COMPLETE FOR EVERY SEQUENCING UNIT, i.e., lane for Illumina data, chip for Ion Torrent data

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "Run this in the desired output directory.  All realigned BAM files in the specified directory (below) will be analysed."
printf '\n\n'

echo "What directory would you like to run GATK Base Recalibration on? This should be done on a per-chip or per-lane basis."
#read -e rootdir
rootdir=/home/mmeier/DNA/data/mutamouse/variant_calling/alignments/raw/bwa_7.12/
echo $rootdir
echo "What directory would you like to put the output files in? Please provide the full path."
#read -e outdir
outdir=/home/mmeier/DNA/data/mutamouse/variant_calling/alignments/realigned
echo $outdir
echo "Would you like to grep for a particular pattern?"
read pattern

files=$(find $rootdir -name "*_dedup.bam" -type f | grep -v old | sort)
processed=$(find $outdir -name "*.recal_reads.bam" -type f | sort)

# Count number of files in input; report them as raw filenames
num=$(for i in $files; do echo $i; done | wc -l)
echo "Selected input directory has the following $num files:"
# Create an array of input files
infileArray=()
for i in $files
do
echo $i
# Substitute for base filename
inputFile=${i##*/}
inputFile=${inputFile%%_dedup.bam}
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
doneFile=${doneFile%%.recal_reads.bam*}
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
#runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) | sed "s|^|$rootdir\/|" | sed "s|$|\.bam.realigned_reads.bam|" )

# Cheats to add root directory and extension: sed "s|^|$rootdir\/|" | sed "s|$|\.bam|"
# THE CODE ABOVE COULD BE REPLACED WITH A WAY TO ASSIGN THE PROPER DIRECTORY STRUCTURE - This would take some work though.


runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) )
# find $rootdir -name "*.realigned_reads.bam" -type f | grep -v old | grep -v analysis_files | sort

FILEPATHS=()
for i in $runlist
do
j=$(find $rootdir -name "*$i*_dedup.bam" -type f | grep -v old )
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
fname=${fname%%.realigned_reads.bam*}
echo "Starting job for $fname..."

screen -d -m -S "$fname".base_recal bash -c "

java -jar $GATK36 \
    -T BaseRecalibrator \
    -R $mm10	 \
    -I $i \
    --knownSites $mm10_snps \
    --knownSites $mm10_indels \
    -o $outdir/$fname.recal_data.table;

java -jar $GATK36 \
    -T BaseRecalibrator \
    -R $mm10 \
    -I $i \
    --knownSites $mm10_snps \
    --knownSites $mm10_indels \
    -BQSR $outdir/$fname.recal_data.table \
    -o $outdir/$fname.post_recal_data.table;

java -jar $GATK36 \
    -T AnalyzeCovariates \
    -R $mm10 \
    -before $outdir/$fname.recal_data.table \
    -after $outdir/$fname.post_recal_data.table \
    -plots $outdir/$fname.recalibration_plots.pdf;

java -jar $GATK36 \
    -T PrintReads \
    -R $mm10 \
    -I $i \
    -BQSR $outdir/$fname.recal_data.table \
    -o $outdir/$fname.recal_reads.bam 
"


# Make sure things start ok
sleep 10

# Don't start more than X instances of screen at once
numProc=$(ps auxw|grep -i GenomeAnalysis|grep -v grep|wc -l)
while [ $numProc -gt 40 ]
do sleep 5
numProc=$(ps auxw|grep -i GenomeAnalysis|grep -v grep|wc -l)
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
