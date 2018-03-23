#!bin/sh
# This script will loop through all the sequences from a BAM header from an input file, and run a separate thread for each sequence.  This is useful for running samtools in a somewhat multithreaded manner.

echo "What directory would you like to run samtools mpileup on?"
read -e rootdir

bamFiles=$(find $rootdir -name "*.bam")

echo "The BAM files to be processed are:"
echo $bamFiles

echo "The BAM files, line by line, are:"
for i in $bamFiles
do
echo $i
bamArray+=($i)
done

# Store header from BAM file
# NOTE they must be aligned to the same reference sequence, because we are only taking the header from ONE
echo "Grabbing header from ${bamArray[0]}"
bamHeader=$(samtools view -H ${bamArray[0]})

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
	# You have 1 second to kill it
	sleep 1

###########################################################
###########################################################
###  START COMMANDS #######################################
###########################################################
###########################################################

# Write a file containing names of bam files to be processed
printf "%s\n" "${bamArray[@]}" > bamFiles.txt

# Loop through each sequence ID
for id in $myList
do
echo "Starting samtools mpileup for $id..."
echo "Chromosome is $id"
echo "The BAM files to be run are $bamFiles"

# Do a pileup on each chromosome for all the BAM files
screen -d -m -S "$id".samtools bash -c "/utils/appz/samtools/samtools-1.2/samtools mpileup -g -t DP -f $mm10min -b bamFiles.txt -r $id > $id.bcf"
# This step may be better done using merged BAM files for each sample (still split by chromosome)

# nohup samtools mpileup -g -t DP -f /home/manager/data/mm10min.fa $bamFiles -r $id > $id.output.bcf &
# If you use nohup, change this above code to make sure IO redirection is fixed
# Need to redirect stdout and stderr of samtools to a different file/log - otherwise the BCF format will be corrupted by the plaintext messages
# The -o flag could be used to write the output file instead

done

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

