#!/bin/bash

# Takes a directory of VCF files (each containing at least 3 biological samples)
# Runs BCF tools to call SNPS
# Runs GATK PhaseByTransmission to determine family characteristics

printf '\n\n'
echo "PLEASE NOTE THIS SCRIPT MUST BE RUN USING BASH.  NOT SH.  IF YOU USED SH, START AGAIN."
printf '\n\n'
echo "This script will take a directory of VCF files (PRODUCED BY MPILEUP), preferably split by chromosome with each containing SNP information for at least one trio.  It will call variants using BCF Tools and run GATK PhaseByTransmission, optionally outputting a mendelian violations file to look for de novo mutations."
printf '\n\n'
echo "REMOVE -maxdepth 1 ARGUMENT FROM FIND TO DESCEND INTO DIRECTORIES, IF DESIRED.  Currently, the check for processed files parses filenames to the first period; ensure there are no abmiguous names in your directory before proceeding."

echo "What directory would you like to process? IT SHOULD CONTAIN VCFS PRODUCED BY MPILEUP."
read -e rootdir
echo "What directory would you like to put the output files in?"
read -e outdir
echo "What sample would you like to process?  This will grep for the pattern in the filename."
read pattern


files=$(find $rootdir -maxdepth 1 -name "*.bcf" -type f | grep $pattern | sort)


for i in $files
do

fname=${i##*/}
fname=${fname%%.bcf}

echo $i
echo $fname

screen -d -m -S "$fname".BCFtools bash -c "
$bcftoolsexp call -g -m $i > $outdir/$fname.raw
"

# Don't start more than X instances of screen at once
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
while [ $numProc -gt 35 ]
do
sleep 5
numProc=$(ps auxw|grep -i SCREEN|grep -v grep|wc -l)
done

done
