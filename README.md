# FrustraEvo: A Tool To Study Energetic Patterns Within And Between Protein Families

FrustraEvo is a very simple tool to calculate the frustration logo and the frustration contact maps. As a input files receives a set of aligned sequences (MSA) and their corresponding structures and calculates the degree of frustration conservation across them, per position or contact using the frustratometeR.

Dependencies:

hmmer (`sudo apt-get install -y hmmer`)

frustratometeR (https://github.com/proteinphysiologylab/frustratometeR)

R library : argparse, ggplot2, ggseqlogo, cowplot, seqinr and data.table

# Minimum code to calculate frustration in a protein (You need to create a .py file (e.g run_logo.py) and put the following code)
```

import os

import sys

sys.path.append('')#Path to Functions.py file

import Functions

import argparse

parser = argparse.ArgumentParser(description='Calculation of the frustration logo.')

parser.add_argument("--JobId", help="The name of the job")

parser.add_argument("--RPath", default='Scripts', help="Path to R script files (Default: Scripts)")

parser.add_argument("--fasta", help="Name of the fasta")

parser.add_argument("--ref", default='None', help="Id of the reference protein of your logo (Default: None)")

parser.add_argument("--pdb_db", default='None', help="Path to the PDBs folder (Default: None)")

parser.add_argument("--cmaps", default='None', help="Put 'yes' for contactmaps calculation (Default: None)")

args = parser.parse_args()

list_file=''

print('Preparing Files...')

list_file=Functions.PreparingFiles(args.JobId,args.RPath,args.pdb_db,args.fasta)

print('Running Frustration')

Functions.FrustraPDB(list_file,args.JobId,args.pdb_db)

print('Making the plots')

Functions.plots_logo(args.JobId,args.ref,args.RPath)

```

To run in terminal: `python3 run_logo.py --JobId XXX --fasta XXX.fasta --ref XXXXX --pdb_db XXX`

**Sample Input files are available at:**
https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples_InputData

**Sample Output files are available at:**
https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples_OutputData

**Full input and output files to reproduce results from the Freiberger-Ruiz Serrra et al 2023 article (https://www.biorxiv.org/content/10.1101/2023.01.25.525527v1) are available at:**
https://github.com/proteinphysiologylab/Freiberger-RuizSerra_et_al2023

You can find the description of the FrustraEvo results outputs in this README file
https://github.com/proteinphysiologylab/FrustraEvo/blob/master/wiki/README

Enjoy :)
