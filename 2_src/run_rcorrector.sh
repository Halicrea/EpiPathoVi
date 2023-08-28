#!/usr/bin/env bash
######################################
#		05/07/2023
#		By Elyna Bouchereau
######################################

export reads_path="RAW_DATA/EpiVib-RNAseq_Vibrio-aestu12-016_Sel-Cu_Agnes-TRAVERS/"
run_rcorrector.pl -t 40 -k 25 -s \
	${reads_path}15-1_S3_L001_R1_001.fastq.gz,${reads_path}15-2_S6_L001_R1_001.fastq.gz,${reads_path}15-3_S10_L001_R1_001.fastq.gz,${reads_path}22-1_S14_L001_R1_001.fastq.gz,${reads_path}22-2_S17_L001_R1_001.fastq.gz,${reads_path}22-3_S21_L001_R1_001.fastq.gz,${reads_path}34-1_S4_L001_R1_001.fastq.gz,${reads_path}34-2_S7_L001_R1_001.fastq.gz,${reads_path}34-3_S11_L001_R1_001.fastq.gz,${reads_path}40-1_S15_L001_R1_001.fastq.gz,${reads_path}40-2_S18_L001_R1_001.fastq.gz,${reads_path}40-3_S22_L001_R1_001.fastq.gz,${reads_path}XC1_S2_L001_R1_001.fastq.gz,${reads_path}XC2_S5_L001_R1_001.fastq.gz,${reads_path}XC3_S9_L001_R1_001.fastq.gz,${reads_path}XM1_S12_L001_R1_001.fastq.gz,${reads_path}XM2_S16_L001_R1_001.fastq.gz,${reads_path}XM3_S19_L001_R1_001.fastq.gz,${reads_path}XS1_S13_L001_R1_001.fastq.gz,${reads_path}XS2_S24_L001_R1_001.fastq.gz,${reads_path}XS3_S20_L001_R1_001.fastq.gz,${reads_path}XT1_S1_L001_R1_001.fastq.gz,${reads_path}XT2_S23_L001_R1_001.fastq.gz,${reads_path}XT3_S8_L001_R1_001.fastq.gz -od RCORRECT_OUT/ -verbose
