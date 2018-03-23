#!/bin/bash

# Takes a directory full of BAM files
# Runs GATK Realign Target Creator
# COMPLETE FOR EVERY SEQUENCING UNIT, i.e., lane for Illumina data, chip for Ion Torrent data
printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "Run this in the desired output directory.  All BAM files in the specified directory (below) will be analysed."
printf '\n\n'

echo "What directory would you like to run GATK Realigner on? This should be done on a per-chip or per-lane basis."
printf '\n'
echo "Note: subdirectories will also be searched unless this is disabled in the find command by setting maxdepth to 1"
read -e rootdir
echo "What directory would you like to put the output files in? Please provide the full path."
read -e outdir

files=$(find $rootdir -name "*.bam" -type f | grep -v old | grep -v analysis_files | sort)
processed=$(find $outdir -name "*.bam*" -type f | sort)

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
runlist=()
runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) )

# THE CODE ABOVE COULD BE REPLACED WITH A WAY TO ASSIGN THE PROPER DIRECTORY STRUCTURE, SUCH AS THE FOLLOWING:
# Cheats to add root directory and extension: sed "s|^|$rootdir\/|" | sed "s|$|\.bam|"
# Instead, I am using the basename to search for the original files.

FILEPATHS=()
for i in $runlist
do
j=$(find $rootdir -name "*$i*.bam" -type f | grep -v old | grep -v analysis_files)
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
fname=${fname%%.*}
echo "Now starting screen for sample $fname, file path is $i"
echo "The input directory is $rootdir"
echo "The output directory is $outdir"

# Start a screen session detatched and named by the file being processed
# Run commands in new bash shell
screen -d -m -S "$fname".$chr.realign bash -c "
java -jar $GATK \
	-T RealignerTargetCreator \
	-R $mm10min \
	-I $i \
	-nt 1 \
	-o $outdir/$fname.target_intervals.list \
	-log $outdir/$fname.log;
	
java -jar $GATK \
	-T IndelRealigner \
	-R $mm10min \
	-I $i \
	-targetIntervals $outdir/$fname.target_intervals.list \
	-o $outdir/$fname.realigned_reads.bam \
	-log $outdir/$fname.log

"

# Could split RTC by chromosome, as well.

# Make sure things start ok
sleep 10

# Don't start more than 30 instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 20 ]
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
