import sys
from Bio import SeqIO

pathAlign=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign4.fasta"
salida1=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SalidaAlign.fasta",'w')
salida2=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SalidaAlignSE.fasta",'w')
salida3=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Posiciones",'w')

vector = list()
c=0	
m=0
control=0 
cgaps=0;
tam=0
tam=0
r=0

for seq_record in SeqIO.parse(pathAlign, "fasta"):
	seqid=seq_record.id
	seq=seq_record.seq
	pdbid=seqid.split("_")
	sres=""
	longline=0
	if c!=0:
		if len(pdbid) >=2:
			pdb=pdbid[0]+"_"+pdbid[1]
		else:
			pdb=pdbid[0]
		salida1.write(">"+pdb+"\n")
		salida2.write(">"+pdb+"\n")
		salida3.write(">"+pdb+"\n")
		sres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+pdb+".done/FrustrationData/"+pdb+".pdb_singleresidue",'r')
		lsres=sres.readline()
	if c==0:
		c=1
		tam=len(seq)
		i=0
		while(i<tam):
			if seq[i] == "-" or seq[i] == "Z":
				vector.append(0)
			else:
				vector.append(1)
			i += 1
	else:
		splitres=""
		tam=len(seq)
		lsres = sres.readline()
		if r==1 and len(pdbid)>2:
			splitres=lsres.split(" ")
			r= int(splitres[0])
			a = int(pdbid[2])-1
			while (r<a):
				lsres=sres.readline()
				r = r + 1
		n=0
		q=0
		for j in range (0,tam):
			if vector[j] == 0:
				if seq[j] != "-":
					lsres = sres.readline()
			else:
				if seq[j] == "Z":
					salida1.write("-")
					longline+=1
					n = n + 1
					q = q + 1
					salida2.write("-")
				else:
					salida1.write(str(seq[j]))
					longline+=1
					n = n + 1
					q = q + 1
					salida2.write(str(seq[j]))
				if n==51:
					salida1.write("\n")
					n=0
				if seq[j] == "Z":
					salida3.write("Z ")
				elif seq[j] == "-":
					salida3.write("G ")
				else:
					if len(splitres) < 2:
						splitres = lsres.split(" ")
					lsres = sres.readline()
					salida3.write(str(splitres[0])+" ")
					splitres = lsres.split(" ")
					r = r + 1
					
		salida1.write("\n")
		salida2.write("\n")
		salida3.write("\n")
sres.close()


longalign=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Equivalences/long.txt",'w')
longalign.write(str(longline))
longalign.close()

salida1.close()
salida2.close()
salida3.close()

