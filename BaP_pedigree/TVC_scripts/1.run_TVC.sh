#!/bin/bash

echo "What is the root directory where your BAM files are stored?"
read -e rootdir
echo "$rootdir"

echo "What directory would you like to run this script in? Note it will create a new directory for each sample being processed."
echo "PLEASE PAY ATTENTION, THE FILES WILL BE WRITTEN HERE:"
read -e outdir
echo "$outdir"

read -r -p "Please ensure the above is correct before you continue!!! [y/N] " response
case $response in
    [yY][eE][sS]|[yY]) 

        echo "Proceeding to next step."
	# You have 2 seconds to kill it
	sleep 2
		

echo "Please provide a filename containing on each line the samples you would like to process"
read -e sampleList

echo "Please provide a file containing JSON formatted parameters to pass to the TVC"
# read -e params
params=germline_low_stringency_proton.json
echo $params

##### CODE BLOCK #####

# SAMPLE LOOP
for i in $(cat $sampleList)

do

# GO TO OUTPUT DIRECTORY, MAKE NEW DIRECTORY FOR EACH SAMPLE AND ENTER IT
cd $outdir
mkdir $i
cd $i

echo "Sample $i lists the following BAM files:"

# Get file list for sample
files=$(find $rootdir -name *_"$i"_*.bam | tr '\n' ',' | rev | cut -c 2- | rev)
bamFile=$(find $rootdir -name *_"$i"_*.bam | head -n 1)

# List files found for sample
for j in $files; do echo $j; done

# Store header from BAM file
# NOTE they must be aligned to the same reference sequence, because we are only taking the header from ONE
echo "Grabbing header from $bamFile... Found the following chromosomes:"
bamHeader=$(samtools view -H $bamFile)

# Obtain a list of all the sequences from the BAM headers
# Grab lines starting with @SQ, then take the second column and used sed to remove the SN: part
# chrList=$(echo "$bamHeader" | awk -F" " '$1 == "@SQ"{print $2}' | sed 's/SN\://g')
# for k in $chrList
# do
# echo "Starting TVC for $i $k"
# touch $i.$k.test

echo "Starting TVC for $i..."
#### TVC RUN CODE ####

# /utils/appz/TVC/tvc-5.0.2-CentOS_6.6_x86_64-binary/bin/variant_caller_pipeline.py -i $files -r $mm10min --num-threads=36 --generate-gvcf=on --region-bed=/home/ngs/references/mouse/ucsc/mm10/index/samtools/non_repeat_regions.bed -p tvc_params.txt
screen -d -m -S "$i".TVC bash -c "
python /utils/appz/TVC/tvc-5.0.2-CentOS_6.6_x86_64-binary/bin/variant_caller_pipeline.py -i $files -r $mm10min --num-threads=14 --generate-gvcf=on -p $params
"

#### END TVC ####


# While running many TVC threads... don't start too many
numProc=$(ps auxw|grep -i variant_caller_pipeline.py|grep -v grep|wc -l)
while [ $numProc -gt 5 ]
do sleep 5
numProc=$(ps auxw|grep -i variant_caller_pipeline.py|grep -v grep|wc -l)
done

# done

done





        ;;
    *)
        echo "Quitting now"
        ;;
esac

