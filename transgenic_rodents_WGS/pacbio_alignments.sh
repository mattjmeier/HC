
# Do alignment of fastq files to mm10 based on input and output folders

# Run line by line.... not a real script!

in=/s0/ngs/DNA/ref_gen/mutamouse/pacbio
out=/home/mmeier/DNA/data/mutamouse/variant_calling/alignments/pacbio

PU=GQPACBIO
SM=MutaMouse14
PL=PACBIO

ID=Run154_CCS
LB=1
files=$(find $in -name "*.ccs.fastq" | grep 154)

ID=Run155_CCS
LB=1
files=$(find $in -name "*.ccs.fastq" | grep 155)

ID=Run246_CCS
LB=1
files=$(find $in -name "*.ccs.fastq" | grep 246)

ID=Run154_SUBREADS
LB=1
files=$(find $in -name "*.filtered_subreads.fastq" | grep 154)

ID=Run155_SUBREADS
LB=1
files=$(find $in -name "*.filtered_subreads.fastq" | grep 155)

ID=Run246_SUBREADS
LB=1
files=$(find $in -name "*.filtered_subreads.fastq" | grep 246)



for i in $files
do echo $i
fname=${i##*/}
fname=${fname%%.fastq}
echo $fname
echo $ID
screen -d -m -S "$ID"."$SM".bwa-mem bash -c "/utils/appz/bwa/bwa.kit/bwa mem -t 1 -R '@RG\tID:"$ID"\tSM:"$SM"\tPL:"$PL"\tLB:"$LB"\tPU:"$PU"' -M -x pacbio $mm10min $i  > $out/$fname.sam"
done






samtools view -Sbh
samtools sort $out/$fname.bam
samtools index $out/$fname.sorted.bam


