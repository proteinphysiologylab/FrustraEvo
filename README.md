# FrustraEvo: A Tool To Study Energetic Patterns Within And Between Protein Families

FrustraEvo is a very simple tool to calculate the frustration logo and the furstration contact maps. As a input files receives a set of aligned sequences (MSA) and their corresponding structures and calculates the degree of frustration conservation across them, per position or contact using the frustratometeR.

Dependencies:

hmmer

frustratometeR (https://github.com/proteinphysiologylab/frustratometeR)

# Minimum code to calculate frustration in a protein (You need to creat a .py file and put the code below)

`import os\n
 import sys
 sys.path.append('')#Path to Functions.py file
 import Functions
 jodib=sys.argv[1]
 path_to_r=sys.argv[2]
 fasta_file=sys.argv[3]
 list_file=''
 prot_ref=sys.argv[4]
 path_to_Pdbs=sys.argv[5]
 Functions.copyfiles(jodib,path_to_r,path_to_Pdbs)
 list_file=Functions.pdb_list(fasta_file)
 Functions.changes(jodib,fasta_file)
 Functions.FrustraPDB(list_file,jodib,path_to_Pdbs)
 Functions.checks(jodib)
 Functions.prepare_file(jodib,prot_ref)
 Functions.FinalAlign(jodib)
 Functions.Equivalences(jodib)
 Functions.LogoCheck(jodib)
 Functions.plots_logo(jodib,prot_ref)`

## **You can find an example of how to use the package at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples

## **You can also find useful examples in our wiki!!:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/wiki

## **You can find the data used in the Wiki examples at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Data
