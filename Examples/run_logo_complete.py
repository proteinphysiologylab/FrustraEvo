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

#How to run the pipeline in linux terminal:
#python3 run_logo.py --JobId XXX --fasta XXX.fasta --ref XXX --pdb_db XXX


args = parser.parse_args()
list_file=''

#the parameter Scripts is the path to the .r and .py files for plots

print('Coping files for Frustration Logo')
Functions.copyfiles(args.JobId,args.RPath,args.pdb_db)
list_file=Functions.pdb_list(args.fasta)
print('Changes in MSA Frustration Logo')
Functions.changes(args.JobId,args.fasta)
print('Running Frutration')
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
print('Making Visualization scripts (.pml)')
Functions.VScript(args.JobId)
print('Running CMaps for Mutational')
Functions.CMaps_Mutational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice mutational
print('Running CMaps for Configurational')
Functions.CMaps_Configurational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice configurational
