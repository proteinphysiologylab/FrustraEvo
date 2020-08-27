import sys
import os

lista=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignSE.fasta",'r')
sal=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign.fasta",'w')
salida=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/ListaPDBC.txt",'w')
errorsal=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/EvoFrustra-log",'w')

while True:
	linea=lista.readline()
	if not linea:
		break
	linea=linea[:-1]
	spl=linea.split("_")
	ta=len(spl);
	c=0
	error=0
	if linea[0] == ">":
		line=linea[1:]
		spline=line.split("_")
		if spline[1] == "":
			frustra=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+spline[0]+".pdb.done/FrustrationData/"+spline[0]+".pdb_msingleresidue")
		else:
			frustra=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+spline[0]+"_"+spline[1]+".pdb.done/FrustrationData/"+spline[0]+"_"+spline[1]+".pdb_msingleresidue")
		for fline in frustra.readlines():
			c=c+1
			spfrus= fline.split(" ")
			tam=len(spfrus)
			c=c+1
			if ta==1:
				if tam<7:
					error=0
				else:
					if c>4:
						error=1
			else:
				if tam<8:
					error=0
				else:
					if c>4:
						error=1
		frustra.close()
		if error==0:
			errorsal.write(linea+'\n')
			linea=lista.readline()
		else:
			line=linea[1:]
			salida.write(line+'\n')
			sal.write(linea+'\n')
			linea=lista.readline()
			sal.write(linea)
lista.close()
salida.close()
sal.close()
errorsal.close()
