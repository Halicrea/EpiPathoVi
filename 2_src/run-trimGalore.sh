#!/usr/bin/env bash
######################################
#		30/06/2023
#		By Elyna Bouchereau
######################################
echo "#### Runing TrimGalore ####"

find SortMeRNA_OUT/ \
	-name "*.fq.gz" | parallel --bar 'Sample_name=$(echo "{}" | cut -f 3 -d"/" | cut -f 1 -d"_"); 
	echo "Working on: ${Sample_name}.";
	trim_galore --illumina -q 30 --phred33 -j 3 --length 35 -o TRIMMED/ --fastqc --fastqc_args "--outdir TRIMMED/FastQC" --three_prime_clip_R1 1 --three_prime_clip_R2 1 "{}" &> ${Sample_name}_trimmed.log'
