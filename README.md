# FrustraEvo: A Tool To Study Energetic Patterns Within And Between Protein Families

FrustraEvo is a very simple tool to calculate the frustration logo and the furstration contact maps. As a input files receives a set of aligned sequences (MSA) and their corresponding structures and calculates the degree of frustration conservation across them, per position or contact using the frustratometeR.

Dependencies:

hmmer

frustratometeR (https://github.com/proteinphysiologylab/frustratometeR)

# Minimum code to calculate frustration in a protein (You need to creat a .py file an put the code below)

`import os`

`import sys`

`sys.path.append('')#Path to Functions.py file`

`import Functions`

`jodib=sys.argv[1]`

`path_to_r=sys.argv[2]`

`fasta_file=sys.argv[3]`

`list_file=''`

`prot_ref=sys.argv[4]`

`path_to_Pdbs=sys.argv[5]`

`print('Coping files for Frustration Logo')`

`Functions.copyfiles(jodib,path_to_r,path_to_Pdbs)`

`list_file=Functions.pdb_list(fasta_file)`

`print('Changes in MSA Frustration Logo')`

`Functions.changes(jodib,fasta_file)`

`print('Running Frutration')`

`Functions.FrustraPDB(list_file,jodib,path_to_Pdbs)`

`print('Running Checks')`

`Functions.checks(jodib)`

`print('Preparing MSA Files to process')`

`Functions.prepare_file(jodib,prot_ref)`

`Functions.FinalAlign(jodib)`

`print('Running Equivalences')`

`Functions.Equivalences(jodib)`

`print('Running Checks')`

`Functions.LogoCheck(jodib)`

`print('Making the plots')`

`Functions.plots_logo(jodib,prot_ref)`

## **You can find an example of how to use the package at:**

https://github.com/proteinphysiologylab/FrustraEvo/Examples

## **You can also find useful examples in our wiki!!:**

https://github.com/proteinphysiologylab/FrustraEvo/wiki

## **You can find the data used in the Wiki examples at:**

https://github.com/proteinphysiologylab/FrustraEvo/Data
