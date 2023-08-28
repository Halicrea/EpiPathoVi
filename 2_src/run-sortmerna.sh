#!/usr/bin/env bash
######################################
#		30/06/2023
#		By Elyna Bouchereau
######################################
echo "#### Runing SortMeRNA ####"
#for samples in RAW_DATA/EpiVib-RNAseq_Vibrio-aestu12-016_Sel-Cu_Agnes-TRAVERS/*.fastq.gz; do
#	Sample_name=$(echo "$samples" | cut -f 3 -d'/' | cut -f 1 -d'_')
#	echo "Working on: ${Sample_name}."
#	sortmerna --ref smr_v4.3_sensitive_db_rfam_seeds.fasta \
#		--reads ${samples} \
#		--workdir SortMeRNA_${Sample_name}_OUT/ \
#		--idx-dir SortmeRNA_DB/	
#done

find RAW_DATA/EpiVib-RNAseq_Vibrio-aestu12-016_Sel-Cu_Agnes-TRAVERS/ -name "*.fastq.gz" | parallel -j 24 'Sample_name=$(echo "{}" | cut -f 3 -d"/" | cut -f 1 -d"_"); echo "Working on: ${Sample_name}.";sortmerna --ref smr_v4.3_default_db.fasta --reads "{}" --workdir SortMeRNA/SortMeRNA_${Sample_name}_OUT/ --idx-dir SortMeRNA/SortMeRNA_DB/ --aligned SortMeRNA/SortMeRNA_${Sample_name}_OUT/${Sample_name}_rRNA-reads --other SortMeRNA/SortMeRNA_${Sample_name}_OUT/${Sample_name}_non-rRNA-reads --fastx --threads 2'
