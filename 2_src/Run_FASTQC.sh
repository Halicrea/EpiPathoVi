#!/usr/bin/env bash
######################################
#		28/06/2023
#		By Elyna Bouchereau
######################################

#Run fastqc on everything in parallel.
find $1 -name '*.fq.gz' | awk '{printf("fastqc \"%s\"\n", $0)}' | parallel -j 24 --verbose
# copies all the fastqc files to directory ./    
find $1 -name '*fastqc.*' | xargs -I '{}' mv '{}' ./TRIMMED_non-rRNA/FASTQC_OUT/
