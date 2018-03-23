#############################################################
#############################################################
# Some useful one-liner code and random bits of information #
#############################################################
#############################################################

####################################
## Rename a sample within a BAM file
####################################

samtools view -h $file.bam | sed 's/Name/NewName/g' | samtools view -Sbh - > $outfile.bam

#########################################################################################################
### Change contig lines in a BAM header using the following (i.e., if you need to add chromosome lengths)
#########################################################################################################

##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>
##contig=<ID=chr6>
##contig=<ID=chr7>
##contig=<ID=chr8>
##contig=<ID=chr9>
##contig=<ID=chr10>
##contig=<ID=chr11>
##contig=<ID=chr12>
##contig=<ID=chr13>
##contig=<ID=chr14>
##contig=<ID=chr15>
##contig=<ID=chr16>
##contig=<ID=chr17>
##contig=<ID=chr18>
##contig=<ID=chr19>
##contig=<ID=chrX>
##contig=<ID=chrY>
##contig=<ID=chrM>

##contig=<ID=chr1,length=195471971>
##contig=<ID=chr2,length=182113224>
##contig=<ID=chr3,length=160039680>
##contig=<ID=chr4,length=156508116>
##contig=<ID=chr5,length=151834684>
##contig=<ID=chr6,length=149736546>
##contig=<ID=chr7,length=145441459>
##contig=<ID=chr8,length=129401213>
##contig=<ID=chr9,length=124595110>
##contig=<ID=chr10,length=130694993>
##contig=<ID=chr11,length=122082543>
##contig=<ID=chr12,length=120129022>
##contig=<ID=chr13,length=120421639>
##contig=<ID=chr14,length=124902244>
##contig=<ID=chr15,length=104043685>
##contig=<ID=chr16,length=98207768>
##contig=<ID=chr17,length=94987271>
##contig=<ID=chr18,length=90702639>
##contig=<ID=chr19,length=61431566>
##contig=<ID=chrX,length=171031299>
##contig=<ID=chrY,length=91744698>
##contig=<ID=chrM,length=16299>

# Use sed to replace all the wrong contig lines:

zcat $i | sed 's/##contig=<ID=chr1>/##contig=<ID=chr1,length=195471971>/' |  sed 's/##contig=<ID=chr2>/##contig=<ID=chr2,length=182113224>/' |  sed 's/##contig=<ID=chr3>/##contig=<ID=chr3,length=160039680>/' |  sed 's/##contig=<ID=chr4>/##contig=<ID=chr4,length=156508116>/' |  sed 's/##contig=<ID=chr5>/##contig=<ID=chr5,length=151834684>/' |  sed 's/##contig=<ID=chr6>/##contig=<ID=chr6,length=149736546>/' |  sed 's/##contig=<ID=chr7>/##contig=<ID=chr7,length=145441459>/' |  sed 's/##contig=<ID=chr8>/##contig=<ID=chr8,length=129401213>/' |  sed 's/##contig=<ID=chr9>/##contig=<ID=chr9,length=124595110>/' |  sed 's/##contig=<ID=chr10>/##contig=<ID=chr10,length=130694993>/' |  sed 's/##contig=<ID=chr11>/##contig=<ID=chr11,length=122082543>/' |  sed 's/##contig=<ID=chr12>/##contig=<ID=chr12,length=120129022>/' |  sed 's/##contig=<ID=chr13>/##contig=<ID=chr13,length=120421639>/' |  sed 's/##contig=<ID=chr14>/##contig=<ID=chr14,length=124902244>/' |  sed 's/##contig=<ID=chr15>/##contig=<ID=chr15,length=104043685>/' |  sed 's/##contig=<ID=chr16>/##contig=<ID=chr16,length=98207768>/' |  sed 's/##contig=<ID=chr17>/##contig=<ID=chr17,length=94987271>/' |  sed 's/##contig=<ID=chr18>/##contig=<ID=chr18,length=90702639>/' |  sed 's/##contig=<ID=chr19>/##contig=<ID=chr19,length=61431566>/' |  sed 's/##contig=<ID=chrX>/##contig=<ID=chrX,length=171031299>/' |  sed 's/##contig=<ID=chrY>/##contig=<ID=chrY,length=91744698>/' |  sed 's/##contig=<ID=chrM>/##contig=<ID=chrM,length=16299>/'  > $i.rehead.vcf

# This could be implemented to use reheader, as below...

for i in $files; do zcat $i | head -n 1000 | grep "^#" | sed 's/##contig=<ID=chr1>/##contig=<ID=chr1,length=195471971>/' |  sed 's/##contig=<ID=chr2>/##contig=<ID=chr2,length=182113224>/' |  sed 's/##contig=<ID=chr3>/##contig=<ID=chr3,length=160039680>/' |  sed 's/##contig=<ID=chr4>/##contig=<ID=chr4,length=156508116>/' |  sed 's/##contig=<ID=chr5>/##contig=<ID=chr5,length=151834684>/' |  sed 's/##contig=<ID=chr6>/##contig=<ID=chr6,length=149736546>/' |  sed 's/##contig=<ID=chr7>/##contig=<ID=chr7,length=145441459>/' |  sed 's/##contig=<ID=chr8>/##contig=<ID=chr8,length=129401213>/' |  sed 's/##contig=<ID=chr9>/##contig=<ID=chr9,length=124595110>/' |  sed 's/##contig=<ID=chr10>/##contig=<ID=chr10,length=130694993>/' |  sed 's/##contig=<ID=chr11>/##contig=<ID=chr11,length=122082543>/' |  sed 's/##contig=<ID=chr12>/##contig=<ID=chr12,length=120129022>/' |  sed 's/##contig=<ID=chr13>/##contig=<ID=chr13,length=120421639>/' |  sed 's/##contig=<ID=chr14>/##contig=<ID=chr14,length=124902244>/' |  sed 's/##contig=<ID=chr15>/##contig=<ID=chr15,length=104043685>/' |  sed 's/##contig=<ID=chr16>/##contig=<ID=chr16,length=98207768>/' |  sed 's/##contig=<ID=chr17>/##contig=<ID=chr17,length=94987271>/' |  sed 's/##contig=<ID=chr18>/##contig=<ID=chr18,length=90702639>/' |  sed 's/##contig=<ID=chr19>/##contig=<ID=chr19,length=61431566>/' |  sed 's/##contig=<ID=chrX>/##contig=<ID=chrX,length=171031299>/' |  sed 's/##contig=<ID=chrY>/##contig=<ID=chrY,length=91744698>/' |  sed 's/##contig=<ID=chrM>/##contig=<ID=chrM,length=16299>/'  > $i.new_header; done

# Or use bcftools to rehead the vcf

bcftools reheader -h header.txt -o $kid.trio.merged.reheader.vcf.gz $kid.trio.merged.vcf.gz

java -jar picard.jar SortVcf I=$kid.trio.merged.reheader.vcf.gz O=$kid.trio.merged.sorted.vcf.gz


#################################################################################################################
# To get a list of files with multiple matching arguments - replace the Name1, Name2, Name3 in parentheses below.
#################################################################################################################

find foldername -maxdepth 1 -regextype posix-egrep -regex ".*(Name1|Name2|Name3)_.*.bam" -type f


################################################################################################
# Search for common elements in two lists; see options for comm to output different conditions.
################################################################################################

comm -1 -2 <(cat $file1 | sort ) <(cat $file2 | sort)


##############################
# Change version of VCF file
##############################

cat input.vcf | sed 's/VCFv4.2/VCFv4.1/' | sed 's/,Version=3>/>/' | sed 's/,Version=\"3\">/>/' | sed 's/Number=R/Number=./' > output.vcf


#############################################################################################################################################################################################
# Make output from Variant Effect Predictor (VEP) more concise by removing lines that share common fields, such as transcript ID (this can be changed by altering the columns parsed by awk)
#############################################################################################################################################################################################

for i in *.txt; do echo $i; name=$(echo $i | sed 's/.txt//') ; cat $i |awk '$4==last{next} {last=$4} 1' |  awk '$1==last{next} {last=$1} 1' | awk '$5==last{next} {last=$5} 1' > ./concise/${name}.concise.txt; done

############################################################################################################################################
# Manually create a VCF file from the output of lacZ scripts (i.e., "common calls" file). This is highly inadvisable. There are better ways.
############################################################################################################################################

# cat file.txt | awk ' {print "lacZ\t"$7"\t\.\t"$8"\t"$9"\t1\tPASS\tBG="$12"\tGT:FQ:CCT\t1/1:"$13":"$16} ' ## No header included here. See below.

for i in *.txt; do echo $i; sm=${i%.*}; cat header.txt | sed "s/SAMPLE/$sm/" > $i.vcf ; cat $i | awk ' {print "lacZ\t"$7"\t\.\t"$8"\t"$9"\t1\tPASS\tBG="$12"\tGT:FQ:CCT\t1/1:"$13":"$16} ' >> $i.vcf ; done

##############################################################################################################
# Awk code to split a file by column. Need to modify column numbers and array length depending on your input.
##############################################################################################################

awk '$0 !~ /\#/ {x=split($15,y,"-"); print > y[1]"-"y[2]"-"y[3]}' output.txt


