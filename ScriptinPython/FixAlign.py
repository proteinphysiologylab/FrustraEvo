import sys
from Bio import SeqIO

Lista=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/ListaPDB.txt",'r')
Align=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Alignment.fasta"
Salida=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignSE.fasta",'w')

for seq_record in SeqIO.parse(Align, "fasta"):
	seq=seq_record.seq
	pdbid=Lista.readline()	
	Salida.write(">"+str(pdbid)+str(seq)+"\n")

Lista.close()
Salida.close()
