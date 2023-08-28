#!/usr/bin/env bash
#PBS -q ftp
#PBS -l mem=32g
#PDB -l ncpus=2
#PBS -l walltime=180:00:00
######################################
#               26/07/2023
#               By Elyna Bouchereau
######################################
echo "#######* DOWNLOAD PACBIO FILES *#######"
#######*        LOAD MODULES            *#####

#
# This script illustrates how to get data files using wget.
#

# Put here the directory where to save downloaded file(s)
DATA_DIR=$DATAWORK/EpiVib-Sequencage-PacBIO-EMseq_V-aestuarianus_MARS/

# Enter directory where to save files
cd $DATA_DIR/TMP/

# use wget to get data files.
wget --no-verbose --continue -i $DATA_DIR/list-of-urls.txt

# add more wget commands here for each file to retrieve...ftp://
echo "### DONE ----"
~                                  