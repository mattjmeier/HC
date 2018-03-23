#!/bin/bash

# Takes a directory full of FASTQ files; runs bwa-mem on them according to GATK best practices. 
# COMPLETE FOR EVERY SEQUENCING UNIT, i.e., lane for Illumina data, chip for Ion Torrent data

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "Run this in the desired output directory.  All BAM files in the specified directory (below) will be realigned using GATK."
printf '\n\n'

echo "What directory would you like to process?"
read -e rootdir
echo "What directory would you like to put the output files in? Please provide the full path."
read -e outdir
echo "Would you like to grep for a particular pattern? If so, enter it; if not, leave blank and press enter"
#read pattern

files=$(find $rootdir -name "*.bam" -type f ) # | grep -v old | grep -v analysis_files | grep $pattern | sort)
processed=$(find $outdir -maxdepth 1 -name "*.bam*" -type f | sort)

# Some filename parsing using sed (courtesy of JH) - could be used in creating file arrays
# For input files, remove leading slash
# sed 's/.*\///' |
# For output files, remove leading slash and the file extension
# sed 's/.*\///' | sed 's/\.[^.]*$//' |


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
inputFile=${inputFile%%.bam}
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
doneFile=${doneFile%%.bam}
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

# echo "Output of comm"...
# comm -23 <(printf '%s\n' "${infileArray[@]}") <(printf '%s\n' "${doneArray[@]}")
# runlist=$(comm -23 <(printf '%s\n' "${files[@]}") <(printf '%s\n' "${doneArray[@]}"))

# Takes the files unique to the list of input files
# Ignores files that are present in both lists
# Prepends the input directory and appends the extension to make sure the runlist points to the correct files
runlist=()
runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) )

# Cheats to add root directory and extension: sed "s|^|$rootdir\/|" | sed "s|$|\.bam|"
# THE CODE ABOVE COULD BE REPLACED WITH A WAY TO ASSIGN THE PROPER DIRECTORY STRUCTURE - This would take some work though - below is my current solution.


FILEPATHS=()
for i in $runlist
do
echo "Runlist: $i"
j=$(find $rootdir -name "*$i*.bam" -type f )
FILEPATHS+=($j)
done

num=$(for i in "${FILEPATHS[@]}"; do echo $i; done | wc -l)
echo "The final FILEPATHS list has the following $num files:"
for i in "${FILEPATHS[@]}"
do
echo $i
done



# Store header from BAM file
# NOTE they must be aligned to the same reference sequence, because we are only taking the header from ONE
echo "Grabbing header from ${FILEPATHS[0]}"
bamHeader=$(samtools view -H ${FILEPATHS[0]})

# Obtain a list of all the sequences from the BAM headers
# Grab lines starting with @SQ, then take the second column and used sed to remove the SN: part
chrList=$(echo "$bamHeader" | awk -F" " '$1 == "@SQ"{print $2}' | sed 's/SN\://g')

echo "The sequence IDs for each chromosome are:"
for chr in $chrList
do
echo $chr
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
fname=${fname%%.bam*}
echo "Now starting screen for sample $fname, file path is $i"
echo "The input directory is $rootdir"
echo "The output directory is $outdir"



# Start a screen session detatched and named by the file being processed
# Run commands in new bash shell
screen -d -m -S "$fname".realign bash -c "
java -jar $GATK3446 \
	-T RealignerTargetCreator \
	-R $mm10min \
	-I $i \
	-nt 1 \
	-o $outdir/$fname.target_intervals.list \
	-log $outdir/$fname.log;
	
java -jar $GATK3446 \
	-T IndelRealigner \
	-R $mm10min \
	-I $i \
	-targetIntervals $outdir/$fname.target_intervals.list \
	-o $outdir/$fname.realigned_reads.bam \
	-log $outdir/$fname.log
"

sleep 2 # Make sure things start ok

# Don't start more than 10 instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 10 ]
do sleep 5
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
done # Finish numproc loop

done # Finish file loop

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
