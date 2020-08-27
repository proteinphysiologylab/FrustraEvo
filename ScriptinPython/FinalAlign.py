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

for seq_record in SeqIO.parse(pathAlign, "fasta"):
	seqid=seq_record.id
	seq=seq_record.seq
	tam=len(seq)
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
		sres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+pdb+".pdb.done/FrustrationData/"+pdb+".pdb_msingleresidue",'r')
		lsres=sres.readline()	
	r=1
	if c==0:
		c=1
		i=0
		while(i<tam):
			if seq[i] == "-" or seq[i] == "Z":
				vector.append(0)
			else:
				vector.append(1)
			i = i + 1
	else:
		splitres=""
		if r==1 and len(pdbid)>2:
			lsres = sres.readline()
			splitres=lsres.split(" ")
			r=splitres[0]
			r = int(r)
			a = int(pdbid[2])-1
			while (r<a):
				lsres=sres.readline()
				r = r + 1
		j=0
		n=0
		q=0
		lsres = sres.readline()
		while(j<tam):
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
					longline+=1
					n=0
				if seq[j] == "Z":
					salida3.write("Z ")
					lsres = sres.readline()
				else:
					if seq[j] == "-":
						salida3.write("G ")
					else:
						splitres = lsres.split(" ")
						r = r + 1
						salida3.write(str(splitres[0])+" ")
						lsres = sres.readline()
			j = j + 1
		salida1.write("\n")
		salida2.write("\n")
		salida3.write("\n")
sres.close()


longalign=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Equivalences/long.txt",'w')
print(str(longline))
longalign.write(str(longline))
longalign.close()

salida1.close()
salida2.close()
salida3.close()

