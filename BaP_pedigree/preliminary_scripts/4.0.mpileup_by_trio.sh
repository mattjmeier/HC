#!bin/sh
# This script will loop through all the sequences from a BAM header from an input file, and run a separate thread for each sequence.  This is useful for running samtools in a somewhat multithreaded manner.
# Sample name must be specified
# Only one sample at a time, currently.  Must keep track manually.

echo "What directory would you like to run samtools mpileup on?"
# read -e rootdir
rootdir="/s0/ngs/DNA/pedigrees/BaP/raw/"
echo "What sample names should be processed? Please provide a list of exact names, separated by the pipe symbol, |"
read samples
echo "What should we name this family?"
read family


bamFiles=$(find $rootdir -maxdepth 1 -regextype posix-egrep -regex ".*($samples).*.bam")

#echo "What sample would you like to process? Please provide the exact sample name as it would appear in the @RG tag of the BAM file."
#read -e sampleName

# Subset the bamFiles list to include only those with the sample name of interest
echo "All the BAM files in your input directory, line by line, are:"
for i in $bamFiles
do
echo $i
# Grab sample name from BAM file
#sm=$(samtools view -h $i | head -n 100 | sed -n -e 's/^.*SM:\([0-9A-Za-z:]*\).*$/\1/p' | uniq)
# If BAM file sample name matches the requested sample name, add to array
#if [ "$sm" == "$sampleName" ]; then
bamArray+=($i)
#fi
done

echo "The BAM files to be processed for $sampleName are:"
echo ${bamArray[@]}

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
dateTag=$(date +%Y.%m.%d)
printf "%s\n" "${bamArray[@]}" > bamFiles.$family.$dateTag.txt

# Store bam files to be processed in a variable
files=$(printf " -I %s" "${bamArray[@]}")

# Loop through each sequence ID
for id in $myList
do
echo "Starting samtools mpileup for $id..."
echo "The sample ID is $sampleName"
# echo "The BAM files to be run are ${bamArray[@]}"
printf "%s\n" "${bamArray[@]}"

# Do a pileup on each chromosome for all the BAM files
screen -d -m -S "$family"."$id".mpileup bash -c "samtools mpileup -g -t DP -f $mm10min -b bamFiles.$sampleName.$dateTag.txt -r $id > $family.$id.bcf"

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


