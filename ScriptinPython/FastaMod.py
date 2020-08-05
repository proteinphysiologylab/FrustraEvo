import sys
from Bio import SeqIO

fasta=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SalidaAlignSE.fasta"
archisal=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+"Equivalences/Logo.fasta",'w')

for seq_record in SeqIO.parse(fasta, "fasta"):
	seqid=seq_record.id
	seq=seq_record.seq
	archisal.write(str(seq)+"\n")
archisal.close()
