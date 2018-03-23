#!bin/bash
# This script will take a directory and sample name as input
# Searches for BAM files in the diretory
# Runs the Torrent Variant Caller on BAM files


echo "What directory would you like to run GATK Unified Genotyper on?"
read -e rootdir
echo "Where should the output go?"
read -e outdir

bamFiles=$(find $rootdir -maxdepth 1 -name "*.bam")

echo "What sample would you like to process? Please provide the exact sample name as it would appear in the @RG tag of the BAM file."
read -e sampleName

# Subset the bamFiles list to include only those with the sample name of interest
echo "All the BAM files in your input directory, line by line, are:"
for i in $bamFiles
do
echo $i
# Grab sample name from BAM file
sm=$(samtools view -h $i | head -n 100 | sed -n -e 's/^.*SM:\([0-9A-Za-z:]*\).*$/\1/p' | uniq)
# If BAM file sample name matches the requested sample name, add to array
if [ "$sm" == "$sampleName" ]; then
bamArray+=($i)
fi
done

echo "The ${#bamArray[@]} BAM files to be processed for $sampleName are:"
printf "%s\n" "${bamArray[@]}"

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

# Store file names to be processed from infileArray as a list of -I inputs for TVC
files=$(printf " -I %s" "${bamArray[@]}")

# Loop through each sequence ID
for id in $myList
do
echo "Starting GATK Unified Genotyper for $id..."
echo "The sample ID is $sampleName"
# echo "The BAM files to be run are ${bamArray[@]}"
printf "%s\n" "${bamArray[@]}"







# Don't start more than X instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 35 ]
do sleep 1
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
done

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


