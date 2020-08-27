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
	if sseqid >=2:	
		archisres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sseqid[0]+"_"+sseqid[1]+".pdb.done/FrustrationData/"+sseqid[0]+"_"+sseqid[1]+".pdb_msingleresidue",'r')
	else:
		archisres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sseqid[0]+".pdb.done/FrustrationData/"+sseqid[0]+".pdb_msingleresidue",'r')
	c=len(seq)
	d=len(sseqid)
	if(d>2):
		c=sseqid[2]
	archisres.readline()
	l = archisres.readline()
	sl = l.split(" ")
	i=sl[0]
	while(i<c):
		l = archisres.readline()
		i = i + 1
	w=0
	while True:
		srespos =  l.split(" ") 
		if not srespos:
			break
		if seq[w] == '-' or seq[w] == 'X':
			archisal3.write("-")
			archisal4.write("-")
		else:
			if srespos[4] == "Missing":
				archisal3.write("-")
				archisal4.write("Z")
				l = archisres.readline()
			else:
				archisal3.write(seq[w])
				archisal4.write(seq[w])
				l = archisres.readline()				
		w = w + 1
		if w == len(seq):
			break
	archisal4.write("\n")
	archisal3.write("\n")
archisal3.close()
archisal4.close()
archisres.close()
