#/bin/bash

# TAKE A FOLDER OF VCF FILES AS INPUT
# RETURN A USEFUL TABLE WITH WHATEVER PARAMETERS YOU LIKE

files=$(find $PWD -name "*.vcf" | grep filtered)

ref="/home/ngs/references/lacz/lacZFlncCorr/ReferenceGenomeLacZ.fasta"

for vcf_file in $files

do

# Grab file name
fname=${vcf_file%.vcf}
fname=${fname##*/}

# Start each process
echo "Starting screen for $vcf_file"

# Open a screen and run the command
screen -d -m -S "$fname". bash -c "java -Xmx1G -jar $GATK -R $ref -T VariantsToTable -V $vcf_file -F POS -GF GT -F QUAL -F AF -F DP -F RO -F AO -F QD -F TYPE -F LEN -F FXX -F HRUN -o $fname.txt --splitMultiAllelic --sample_rename_mapping_file rename.tsv; exit"
sleep 1

done
