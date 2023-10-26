# FrustraEvo: A Tool To Study Energetic Patterns in Protein Families

FrustraEvo is a very simple tool to calculate the frustration logo and the frustration contact maps. As a input files receives a set of aligned sequences (MSA) and their corresponding structures and calculates the degree of frustration conservation across them, per position or contact using the frustratometeR.

A docker containing to run FrustrEvo in any platform without having to deal with dependencies can be accesed from here:
https://hub.docker.com/r/proteinphysiologylab/frustraevo

Dependencies:

hmmer (`sudo apt-get install -y hmmer`)

frustratometeR (https://github.com/proteinphysiologylab/frustratometeR)

R library : argparse, ggplot2, ggseqlogo, cowplot, seqinr and data.table

# Colab Notebook
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/engelberger/FrustraEvo/blob/master/FrustraEvo_Colab.ipynb)

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


#How to run the pipeline in linux terminal:
#python run_logo.py JobId RPath fasta list_file ref pdb_db 

args = parser.parse_args()
list_file=''

#the parameter Scripts is the path to the .r and .py files for plots

print('Copyng files for Frustration Logo')
Functions.copyfiles(args.JobId,args.RPath,args.pdb_db)
list_file=Functions.pdb_list(args.fasta)
print('Changes in MSA Frustration Logo')
Functions.changes(args.JobId,args.fasta)
print('Running Checks in sequence')
Functions.checks_seq(list_file,args.JobId,args.pdb_db)
print('Running Frustration')
Functions.FrustraPDB(list_file,args.JobId,args.pdb_db)
print('Running Checks')
Functions.checks(args.JobId)
print('Preparing MSA Files to process')
Functions.prepare_file(args.JobId,args.ref)
Functions.FinalAlign(args.JobId)
print('Running Equivalences')
Functions.Equivalences(args.JobId)
print('Preparing file for sequences for Sequence logo')
Functions.FastaMod(args.JobId)
print('Running Checks')
Functions.LogoCheck(args.JobId)
print('Making the plots')
Functions.plots_logo(args.JobId,args.ref,args.RPath)
Functions.clean_files(args.JobId,args.RPath,args.ref)
print('Making Visualization scripts (.pml)')
Functions.VScript(args.JobId)

```

To run in terminal: `python3 run_logo.py --JobId XXX --fasta XXX.fasta --ref XXXXX --pdb_db XXX`

Example with Sample data once you are located inside of the folder that contains both the MSA (alphas.fasta) and the Pdbs containing folder (Pdbs): 

`python3 run_logo.py --JobId 2022122114524431133 --fasta alphas.fasta --ref 2dn1-A --pdb_db Pdbs`

**Sample Input files are available at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples_InputData

**Sample Output files are available at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples_OutputData

**Full input and output files to reproduce results from the Freiberger-Ruiz Serrra et al 2023 article (https://www.biorxiv.org/content/10.1101/2023.01.25.525527v1) are available at:**

https://github.com/proteinphysiologylab/Freiberger-RuizSerra_et_al2023

You can find the description of the FrustraEvo results outputs in this README file
https://github.com/proteinphysiologylab/FrustraEvo/blob/master/wiki/README

Enjoy :)
