#!/usr/bin/env bash
######################################
#		31/05/2023
#		By Elyna Bouchereau
#       Mise en place du pipeline SMRT
######################################

##*** ENVIRONNEMENT DE TRAVAIL ***##
## Creation of the conda environment
conda create -n smrt-env python=3.9
conda activate smrt-env

## SMRT-tools
conda install -c hcc smrtlink-tools
#  Prepare bam file for genome assembly
conda install pbtk
parallel --bar pbindex {} ::: *.bam
parallel --bar bam2fastq -o FASTQ/{} {} ::: *.bam

#  Genome assembly (https://github.com/PacificBiosciences/pbipa)
conda install pbipa

#  ipdSummary

#  pbvalidate

#  summarizeModifications (generates a GFF file)
