##----Cambios en el Alineamiento---
import sys
from Bio import SeqIO

pathAlign=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+"SeqAlign2.fasta"

archisal3=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+"SeqAlign3.fasta",'w')
archisal4=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+"SeqAlign4.fasta",'w')

for seq_record in SeqIO.parse(pathAlign, "fasta"):
	seqid=seq_record.id
	seq=seq_record.seq
	sseqid=seqid.split("_")
	pdbid=sseqid[0]+"_"+sseqid[1]
	archisal3.write(">"+seqid+"\n")	
	archisal4.write(">"+seqid+"\n")		
	archisres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sseqid[0]+"_"+sseqid[1]+".pdb.done/FrustrationData/"+sseqid[0]+".pdb_singleresidue",'r')
	c=len(seq)
	d=len(sseqid)
	if(d>2):
		c=sseqid[2]
	archisres.readline()
	l = archisres.readline()
	sl = l.split(" ")
	i=sl[0]
	print (sl[0])
	while(i<c):
		l = archisres.readline()
		i = i + 1
	p=len(seq)
	w=0
	while(w<p):
		srespos =  l.split(" ")
		if seq[w] == '-' or seq[w] == 'X':
			archisal3.write("-")
			archisal4.write("-")
		else:
			if srespos[1] == "Missing":
				archisal3.write("-")
				archisal4.write("Z")
				l = archisres.readline()
			else:
				archisal3.write(seq[w])
				archisal4.write(seq[w])
				l = archisres.readline()				
		w = w + 1
	archisal4.write("\n")
	archisal3.write("\n")
archisal3.close()
archisal4.close()
archisres.close()
