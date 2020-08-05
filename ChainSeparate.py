##--Separa cadenas en python--##
import sys

pdb=open(sys.argv[3]+"/OutPutFiles"+sys.argv[4]+"/"+sys.argv[1]+".pdb",'r')
archisal=open(sys.argv[3]+"/OutPutFiles"+sys.argv[4]+"/Modeller/"+sys.argv[1]+".pdb",'w')

for linea in pdb.readlines():
	s=linea.split()
	t=len(linea)
	if t > 20:
		if s[0] == "ATOM":
			if linea[21] == sys.argv[2]:
				archisal.write(linea)
		if s[0] == "TER":
			if linea[21] == sys.argv[2]:
				archisal.write(linea)
		if s[0] == "ENDMDL":
			break

pdb.close()
archisal.close()
