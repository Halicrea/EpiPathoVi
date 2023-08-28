#!/usr/bin/env python3
#########################################
## 23/03/2023
## Par Elyna Bouchereau
## Fichier: search_motif.py
###########################################
from typing import Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio import AlignIO
from Bio.Align import AlignInfo
import Bio
from itertools import product
from Bio.Data.IUPACData import ambiguous_dna_values
import sys
#******************************************
input_file = sys.argv[1]

#*************  TEST       ****************
import sys
import os
import re
import argparse
import textwrap
import pprint
from os.path import dirname, basename, abspath, join
#from paramparser import ParamParser

DEBUG = False

IUPAC_DNA = {
             'N': '[ATGC]',
             'W': '[AT]',
             'S': '[CG]',
             'M': '[AC]',
             'K': '[GT]',
             'R': '[AG]',
             'Y': '[CT]',
             'B': '[CGT]',
             'D': '[AGT]',
             'H': '[ACT]',
             'V': '[ACG]',
             '/': '|'
}

DNA_PROB = {
            'N': 1,
            'W': 'DNA_PROB["A"] * 2',
            'S': 'DNA_PROB["G"] * 2',
            'M': 'DNA_PROB["A"] + DNA_PROB["G"]',
            'K': 'DNA_PROB["A"] + DNA_PROB["G"]',
            'R': 'DNA_PROB["A"] + DNA_PROB["G"]',
            'Y': 'DNA_PROB["A"] + DNA_PROB["G"]',
            'B': 'DNA_PROB["A"] + DNA_PROB["G"] * 2',
            'D': 'DNA_PROB["A"] * 2 + DNA_PROB["G"]',
            'H': 'DNA_PROB["A"] * 2 + DNA_PROB["G"]',
            'V': 'DNA_PROB["A"] + DNA_PROB["G"] * 2',
            'A': 0.25,
            'T': 0.25,
            'G': 0.25,
            'C': 0.25
}

#o_param = ParamParser()
#********************************************************

#motif = [Seq("GATC")]

motif = [Seq("GATC"),
    Seq("CCWGG"),
    Seq("TAACNNNNRTAC"),
    Seq("CAGNNNNNNTYTC"),
    Seq("CAGCTG"),
	Seq("ACCNNNNNNNTTCY")]


def ambiguous_dna_list(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))

motifs_list = []
for i in range(len(motif)):
	motifs_list.append(ambiguous_dna_list(motif[i]))
m = motifs.create(ambiguous_dna_list(motif[5]))
len(m)


if len(sys.argv) < 3:
	output_file = "Test.txt"
else:
	output_file = sys.argv[2]

fasta_seqs = SeqIO.parse(open(input_file),'fasta')
with open(output_file,'w') as out_file:
	for fasta in fasta_seqs:
		name, sequence = fasta.id, str(fasta.seq)
		for pos, seq in m.instances.search(sequence):
			out_file.write("%i %s\n" % (pos, seq))

