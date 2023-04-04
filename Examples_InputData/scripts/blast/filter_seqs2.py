#!/usr/bin/python
# Import packages
import os
import sys
# Import packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import pairwise2
from Bio.Seq import Seq
from pyfaidx import Fasta, wrap_sequence


file = sys.argv[1]
output_file = file
print(file)
with open(output_file + ".nosynthetic", 'a') as out:
    fasta_sequences = SeqIO.parse(open(file),'fasta')
    for fasta in fasta_sequences:
        name = fasta.name
        desc = fasta.description
        seq=fasta.seq
        if '>' in desc: 
            desc_s = desc.split('>')[0]
        else:
            desc_s = desc
        if 'synthetic' not in desc_s and 'artificial' not in desc_s and 'fragment' not in desc_s and 'low quality' not in desc_s and 'partial' not in desc_s and 'synthetic' not in desc_s: #and 'X' not in seq:
            with open(output_file + ".lenseq", 'a') as outf:
                outf.write(name + '\t' + str(len(seq)) + '\n')

            out.write('>{}\n'.format( desc_s))
            for line in wrap_sequence(70, str(fasta.seq)):
                out.write(line)