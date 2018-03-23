#!/bin/bash

# Takes a directory full of FASTQ files; runs bwa-mem on them according to GATK best practices. 
# COMPLETE FOR EVERY SEQUENCING UNIT, i.e., lane for Illumina data, chip for Ion Torrent data

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "Run this in the desired output directory.  All FASTQ files in the specified directory (below) will be aligned to mm10 (minimal genome)."
printf '\n\n'

echo "What directory would you like to process?"
#read -e rootdir
rootdir=/s0/ngs/DNA/pedigrees/BaP/illumina
#rootdir=/s0/ngs/DNA/ref_gen/mutamouse/illumina/HCNextSeq
echo $rootdir

echo "What directory would you like to put the output files in? Please provide the full path."
#read -e outdir
outdir=/home/mmeier/DNA/data/mutamouse/variant_calling/alignments/raw/bwa_7.12/
echo $outdir

#echo "Would you like to grep for a particular pattern? If so, enter it; if not, leave blank and press enter"
#read pattern
# This doesn't work but the idea would be to set variable to . if nothing is entered.
#if [ -z "$pattern" ]
#pattern=.
#fi

echo "What prefix should be used for the libraries? E.g., HC for Health Canada, GQ for Genome Quebec. Applies to the batch of files being processed."
read prefix


files=$(find $rootdir -name "*.fastq*" -type f | sort | grep -v Undetermined)
processed=$(find $outdir -name "*am" -type f | sort)

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
inputFile=${inputFile%%.fastq*}
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
doneFile=${doneFile%.*}
# Add doneFile to array
doneArray+=($doneFile)
done
# Report the processed file array
echo "List of processed files:"
echo ${doneArray[@]}

printf '\n\n'
echo "Output from printf of infileArray:"
printf '\n'
printf '%s\n' "${infileArray[@]}"
printf '\n\n'
echo "Output from printf of doneArray:"
printf '\n'
printf '%s\n' "${doneArray[@]}"
printf '\n'

# Takes the files unique to the list of input files
# Ignores files that are present in both lists
# Prepends the input directory and appends the extension to make sure the runlist points to the correct files
runlist=()
runlist=$(comm -23 <(printf '%s\n' "${infileArray[@]}" |sort ) <(printf '%s\n' "${doneArray[@]}" | sort) )

FILEPATHS=()
for i in $runlist
do
j=$(find $rootdir -name "*$i*.fastq*" -type f | grep -v old | grep -v analysis_files | grep -v R2)
FILEPATHS+=($j)
done

echo "The final FILEPATHS list has the following ${#FILEPATHS[@]} files:"

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
###  THIS IS THE CODE TO RUN BWA MEM ######################
###########################################################
###########################################################

for i in "${FILEPATHS[@]}"

do

# Get file names for R1 and R2 pairs
fname=${i##*/}
fname=${fname%%.fastq*}
fname2=$(echo $fname | sed 's/R1/R2/')
i2=$(echo $i | sed 's/R1/R2/')
echo "Now starting screen for sample $fname, file path is $i and it's pair is $fname2, found at $i2"
echo "The input directory is $rootdir"
echo "The output directory is $outdir"



echo "We must create a read group tag for GATK and other tool compatability. Make sure you understand the meaning of all the terms before using them."
printf '\n'
printf '\n'
echo "FILE: $fname"
echo "PAIR: $fname2"
printf '\n'





echo "What is the ID? This should be unique to this FASTQ file as it will apply to every read aligned here. We will use the flowcell ID for this."
#read ID
# grab ID from run name
#ID=${fname#*.}
#ID=${ID%%.*}
if [[ $i =~ "gz" ]]
then
echo "Detected fastq.gz, using zcat"
ID=$(zcat $i | head -n 1 | awk ' BEGIN {FS=":"}; {gsub(/@/,"",$1)} ; {print $3}')
else echo "Uncompressed fastq detected, using cat"
ID=$(cat $i | head -n 1 | awk ' BEGIN {FS=":"}; {gsub(/@/,"",$1)} ; {print $3}')
fi
echo $ID





echo "What is the biological sample name? E.g., mouse3M"
#read SM
# grab sample name from filename
SM=${fname#*.}
SM=${SM%%_*}
echo $SM

echo "what platform was used for sequencing? ILLUMINA, IONTORRENT, PACBIO, etc."
#read PL
PL=ILLUMINA
echo $PL

echo "What is the library name? For illumina data, use insert size and sample name"
#read LB
# create library information
# if fname matches MP, set MP to reflect
if [[ $fname == *"mate"* ]]; then MP="_mp"; echo "Mate pair library detected, using _mp suffix for library name"; else MP=""; fi
# make LB tag
if [[ -z $prefix ]]
then echo "No prefix set, using default as HC (for Health Canada)"
LB=HC_$SM$MP
echo $LB
else
echo "Prefix is $prefix"
LB="$prefix"_"$SM$MP"
echo $LB
fi

echo "What is the sequencing unit? e.g., lane1 (you may have multiple libraries sequenced on the same lane, and possibly have one library sequenced on multiple lanes)"
# read PU
# Platform Unit
#PU=${fname#*.}
#PU=${fname%%.*}
if [[ $i =~ "gz" ]]
then
echo "Detected fastq.gz, using zcat"
PU=$(zcat $i | head -n 1 | awk ' BEGIN {FS=":"}; {gsub(/@/,"",$1)} ; {print $10}')
else echo "Uncompressed fastq detected, using cat"
PU=$(cat $i | head -n 1 | awk ' BEGIN {FS=":"}; {gsub(/@/,"",$1)} ; {print $10}')
fi
PU=$ID.$PU
echo $PU

echo "How many threads should be devoted to processing this file?"
#read numthreads
numthreads=7
echo $numthreads


echo "Starting now....."
#sleep 2

# Start a screen session detatched and named by the file being processed
# Run commands in new bash shell
screen -d -m -S "$ID"."$SM".bwa bash -c "/utils/appz/bwa/bwa.kit/bwa mem -M -t $numthreads -R '@RG\tID:"$ID"\tSM:"$SM"\tPL:"$PL"\tLB:"$LB"\tPU:"$PU"' $mm10 $i $i2 |  /utils/appz/samtools/samtools-1.3/samtools view -bh > $outdir/$fname.bam"

# Make sure things start ok
sleep 5

# Don't start more than X instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 5 ]
do sleep 5
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
done

done

###########################################################
###########################################################
###  END ##################################################
###########################################################
###########################################################


        ;;
    *)
        echo "Quitting now"
        ;;
esac
