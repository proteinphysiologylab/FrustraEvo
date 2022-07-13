import os
import sys

sys.path.append('')#Path to Functions.py file
import Functions


#How to run the pipeline in linux terminal:
#python run_logo.py jodib path_to_r fasta_file list_file prot_ref path_to_Pdbs 
#python run_logo.py SarsCov Scripts/ SarsCov_ali.fasta Lista_SarsCov NP_000129 PdbsSarsCov

#jodib: name of the job, example: SarsCov
#path_to_r: Path to R script files
#fasta_file: name of the fasta, example: SarsCov_ali.fasta
#prot_ref: Id of the reference protein of your logo, example: NP_000129
#path_to_Pdbs: Path to the PDBs folder


jodib=sys.argv[1]
path_to_r=sys.argv[2]
fasta_file=sys.argv[3]
list_file=''
prot_ref=sys.argv[4]
path_to_Pdbs=sys.argv[5]

#the parameter Scripts is the path to the .r and .py files for plots

print('Coping files for Frustration Logo')
Functions.copyfiles(jodib,path_to_r,path_to_Pdbs)
list_file=Functions.pdb_list(fasta_file)
print('Changes in MSA Frustration Logo')
Functions.changes(jodib,fasta_file)
print('Running Frutration')
Functions.FrustraPDB(list_file,jodib,path_to_Pdbs)
print('Running Checks')
Functions.checks(jodib)
print('Preparing MSA Files to process')
Functions.prepare_file(jodib,prot_ref)
Functions.FinalAlign(jodib)
print('Running Equivalences')
Functions.Equivalences(jodib)
print('Preparing file for sequences for Sequence logo')
Functions.FastaMod(jodib)
print('Running Checks')
Functions.LogoCheck(jodib)
print('Making the plots')
Functions.plots_logo(jodib,prot_ref)
print('Making Visualization scripts (.pml)')
Functions.VScript(jodib)
print('Running CMaps for Mutational')
Functions.CMaps_Mutational(jodib,path_to_r,prot_ref)
print('Running CMaps for Configurational')
Functions.CMaps_Configurational(jodib,path_to_r,prot_ref)
