#!bin/sh
# This script will loop through all the sequences from a BAM header from an input file, and run a separate thread for each sequence.  This is useful for running samtools in a somewhat multithreaded manner.

echo "What directory would you like to run samtools mpileup on?"
read -e rootdir
echo "Where would you like the output?"
read -e outdir

bamFiles=$(find $rootdir -name "*recal_reads.bam"  | sort)

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
#printf "%s\n" "${bamArray[@]}" > bamFiles.txt

# Loop through BAM files
#for i in $bamFiles
#do echo $i

# Loop through each sequence ID
for id in $myList
do
echo "Starting samtools mpileup for $id..."
echo "CHROMOSOME: $id"
echo "BAM FILES: $bamFiles"
#echo "The BAM file to be run is $i"
fname=${i##*/}
fname=${fname%%.bam}
echo "FILE NAME: $fname"
echo "FULL PATH: $outdir/$fname.$id.vcf.gz"

# Do a pileup on each chromosome for all the BAM files
# screen -d -m -S "$id".$i.samtools bash -c "/utils/appz/samtools/samtools-1.3/samtools mpileup -g -t DP -f $mm10min $i -r $id -u $i | /utils/appz/samtools/bcftools-1.3/bcftools call -mv -O z -o $i.$id.vcf.gz -"


nohup bash -c "/utils/appz/samtools/samtools-1.3/samtools mpileup -g -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR  -f $mm10min -r $id -u $bamFiles | /utils/appz/samtools/bcftools-1.3/bcftools call -mv -O z -o $outdir/$fname.$id.vcf.gz " &

# Make sure things start ok
sleep 6

# Don't start more than N instances of screen at once
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)
while [ $numProc -gt 100 ]
do sleep 1
numProc=$(ps auxw|grep -i samtools|grep -v grep|wc -l)
done



done
#done

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

