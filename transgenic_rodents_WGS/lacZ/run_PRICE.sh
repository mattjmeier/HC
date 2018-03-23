#!/bin/bash
# j=1;for i in $files; do printf '%sspf '$i' 1000 4000 \\\n' '-'; let j=$j+1; done

PriceTI \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.G02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.F01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.E01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.D02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.H01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.C02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run154_497.G02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.F02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.G01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.A02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.D01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.E02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.C01_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.B02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run154_497.F02_1.filtered_subreads.longest.fastq 1000 4000 \
	-spf /s0/ngs/DNA/ref_gen/mutamouse/pacbio/Run155_498.H02_1.filtered_subreads.longest.fastq 1000 4000 \
	-mp /s0/ngs/DNA/ref_gen/mutamouse/illumina/HI.2502.MutaMouse14.mate.1.fastq /s0/ngs/DNA/ref_gen/mutamouse/illumina/HI.2502.MutaMouse14.mate.2.fastq 3000 \
	-icf /home/mmeier/DNA/data/mutamouse/lacz/combined_lambda_lac.fasta 1 1 3 \
	-o LacZ_PRICE_assembly_Nov23_lambdagt10lacZ_pacbio_long_illumina_mp.fasta \
	-target 85 5 5 5 \
	-nc 20 \
	-a 30 \
	-log v \
	-logf log_Nov23.txt

