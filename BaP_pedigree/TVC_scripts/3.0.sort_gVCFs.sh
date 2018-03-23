#!/bin/bash
# Sort merged gVCFs and replace header prior to MV discovery

#RUN IN DIRECTORY!!!

echo "Give the sample name to be sorted"
read i

zcat $i.trio.merged.vcf.gz | sed 's/##contig=<ID=chr1>/##contig=<ID=chr1,length=195471971>/' |  sed 's/##contig=<ID=chr2>/##contig=<ID=chr2,length=182113224>/' |  sed 's/##contig=<ID=chr3>/##contig=<ID=chr3,length=160039680>/' |  sed 's/##contig=<ID=chr4>/##contig=<ID=chr4,length=156508116>/' |  sed 's/##contig=<ID=chr5>/##contig=<ID=chr5,length=151834684>/' |  sed 's/##contig=<ID=chr6>/##contig=<ID=chr6,length=149736546>/' |  sed 's/##contig=<ID=chr7>/##contig=<ID=chr7,length=145441459>/' |  sed 's/##contig=<ID=chr8>/##contig=<ID=chr8,length=129401213>/' |  sed 's/##contig=<ID=chr9>/##contig=<ID=chr9,length=124595110>/' |  sed 's/##contig=<ID=chr10>/##contig=<ID=chr10,length=130694993>/' |  sed 's/##contig=<ID=chr11>/##contig=<ID=chr11,length=122082543>/' |  sed 's/##contig=<ID=chr12>/##contig=<ID=chr12,length=120129022>/' |  sed 's/##contig=<ID=chr13>/##contig=<ID=chr13,length=120421639>/' |  sed 's/##contig=<ID=chr14>/##contig=<ID=chr14,length=124902244>/' |  sed 's/##contig=<ID=chr15>/##contig=<ID=chr15,length=104043685>/' |  sed 's/##contig=<ID=chr16>/##contig=<ID=chr16,length=98207768>/' |  sed 's/##contig=<ID=chr17>/##contig=<ID=chr17,length=94987271>/' |  sed 's/##contig=<ID=chr18>/##contig=<ID=chr18,length=90702639>/' |  sed 's/##contig=<ID=chr19>/##contig=<ID=chr19,length=61431566>/' |  sed 's/##contig=<ID=chrX>/##contig=<ID=chrX,length=171031299>/' |  sed 's/##contig=<ID=chrY>/##contig=<ID=chrY,length=91744698>/' |  sed 's/##contig=<ID=chrM>/##contig=<ID=chrM,length=16299>/' | bgzip > $i.trio.merged.rehead.vcf.gz

java -jar /utils/appz/picard/picard/dist/picard.jar SortVcf I=$i.trio.merged.rehead.vcf.gz O=$i.trio.merged.sorted.vcf.gz

