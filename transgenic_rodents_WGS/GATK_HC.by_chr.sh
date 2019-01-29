#!/bin/bash

# Takes a directory full of recalibrated BAM files
# Runs GATK Haplotype Caller in DISCOVERY MODE
# COMPLETE FOR EVERY SAMPLE

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."

echo "What directory would you like to run GATK Variant Discovery on? This should be done on a per-sample basis."
read -e rootdir
echo "What sample would you like to run?"
read sample
files=$(find $rootdir -maxdepth 1 -name "*.bam" -type f | grep -v old | sort -V | grep "$sample" )

echo "What directory would you like to put the output files in? Please provide the full path."
read -e outdir


infileArray=()
for i in $files
do
echo $i
# Add to array
sm=$(samtools view -h $i | head -n 100 | sed -n -e 's/^.*SM:\([0-9A-Za-z_:]*\).*$/\1/p' | uniq)
echo "Sample name: $sm"
infileArray+=($i)
done

for i in ${infileArray[@]}
do
echo $i
done



#fileInput=$(printf " -I %s" "${infileArray[@]}")

read -r -p "Please ensure the above is correct before you continue!!! [y/N] " response
case $response in
    [yY][eE][sS]|[yY])

        echo "Going ahead with analysis..."
        # You have 3 seconds to kill it
        sleep 3

###########################################################
###########################################################
###  START PROGRAM ########################################
###########################################################
###########################################################


for i in ${infileArray[@]}
do


fname=${i##*/}
fname=${fname%.*}

echo "Starting for $fname"

#for chr in $(echo chrM) ### FOR TESTING
for chr in $(cat chr_list.txt)
do
echo $chr

screen -d -m -S "$fname"."$chr".HaplotypeCaller.GVCF bash -c "

java -Xmx20g -jar $GATK36 \
    -T HaplotypeCaller \
    -R $mm10 \
    -I $i \
    -L $chr \
    --emitRefConfidence GVCF \
    -o $outdir/$fname.$chr.g.vcf.gz

"
sleep 30

# Don't start more than N instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 13 ]
do sleep 60
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
# Done checking processes
done

# Done chr loop
done

# Done file loop
done

###########################################################
###########################################################
###  END PROGRAM ##########################################
###########################################################
###########################################################

        ;;
    *)
        echo "Quitting now"
        ;;
esac


