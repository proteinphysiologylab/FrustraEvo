import sys
import os
import os.path as path

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
		if len(spline) == 1:
			pathPDB=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+spline[0]+".done/FrustrationData/"+spline[0]+".pdb_singleresidue"
			if path.exists(pathPDB):
				error=1
		else:
			pathPDB=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+spline[0]+"_"+spline[1]+".done/FrustrationData/"+spline[0]+"_"+spline[1]+".pdb_singleresidue"
			if path.exists(pathPDB):
				error=1
		
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
