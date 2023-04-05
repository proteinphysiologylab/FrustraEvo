These Input Data consist in a Multiple Sequence Alignment (MSA, alphas.fasta) where all the sequences that correspond to the PDB files contained in the Pdbs folder are aligned. 

Each sequence within the MSA has the same header as its corresponding PDB file. When building input files to run through FrustraEvo all sequences headers need to have a corresponding PDB structure with the same exact headers as the PDB file name.

Example for the following two sequences:

>1fsx-A
VLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGAKVAAALTKAVEHLDDLPGALSELSDLHAHKLRVDPVNFKLLSHSLLVTLASHLPSDFTPAVHASLDKFLANVSTVLTSKYR
>3cy5-A
VLSAADKSNIQAAWGKVGGHAADYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGAKVANALTKAVGHLDDLPGALSELSDLHAHKLRVDPVNFKLLSHSLLVTLASHLPNDFTPAVHASLDKFLASVSTVLTSKYR

There need to be two PDB files named as follows within the Pdbs folder:
1fsx_A.pdb
3cy5_A.pdb
