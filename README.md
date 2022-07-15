# FrustraEvo: A Tool To Study Energetic Patterns Within And Between Protein Families

FrustraEvo is a very simple tool to calculate the frustration logo and the furstration contact maps. As a input files receives a set of aligned sequences (MSA) and their corresponding structures and calculates the degree of frustration conservation across them, per position or contact using the frustratometeR.

Dependencies:

hmmer

frustratometeR (https://github.com/proteinphysiologylab/frustratometeR)

# Minimum code to calculate frustration in a protein (You need to creat a .py file and put the code below)

`import os\n`
`import sys\n`
`sys.path.append('')#Path to Functions.py file\n`
`import Functions\n`
`jodib=sys.argv[1]\n`
`path_to_r=sys.argv[2]\n`
`fasta_file=sys.argv[3]\n`
`list_file=''\n`
`prot_ref=sys.argv[4]\n`
`path_to_Pdbs=sys.argv[5]\n`

`Functions.copyfiles(jodib,path_to_r,path_to_Pdbs)\n`
`list_file=Functions.pdb_list(fasta_file)\n`
`Functions.changes(jodib,fasta_file)\n`
`Functions.FrustraPDB(list_file,jodib,path_to_Pdbs)\n`
`Functions.checks(jodib)\n`
`Functions.prepare_file(jodib,prot_ref)\n`
`Functions.FinalAlign(jodib)\n`
`Functions.Equivalences(jodib)\n`
`Functions.LogoCheck(jodib)\n`
`Functions.plots_logo(jodib,prot_ref)\n`

## **You can find an example of how to use the package at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Examples

## **You can also find useful examples in our wiki!!:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/wiki

## **You can find the data used in the Wiki examples at:**

https://github.com/proteinphysiologylab/FrustraEvo/tree/master/Data
