import sys
from Bio import SeqIO

AClean=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignClean.fasta",'w')
pathAlign=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+"SeqAlign.fasta"
vector = list()
c=0

for seq_record in SeqIO.parse(pathAlign, "fasta"):
	c = c + 1
	vector.append(seq_record.seq)
	
i=0
j=0
porc=c*0.6 
cgaps=0;
tam=len(vector[0])
aux=""
vectorgaps = list()
while i<tam:
	j=0
	cgaps=0
	while(j<c):
		aux=vector[j]
		if aux[i] == "-":
			cgaps = cgaps + 1
		j = j + 1 
	if cgaps < porc:
		vectorgaps.append(0)
	else:
		vectorgaps.append(1)
	i = i + 1
tam=len(vectorgaps)
c=0;

for seq_record in SeqIO.parse(pathAlign, "fasta"):
	seqid=seq_record.id
	seq=seq_record.seq
	pr=0
	AClean.write(">"+seqid+"\n")
	i=0
	while(i<tam):
		if vectorgaps[i]==0:
			if pr==60:
				AClean.write("\n")
				pr=0
			AClean.write(seq[i])
			pr+=1
		i = i + 1
	AClean.write("\n")
AClean.close()
