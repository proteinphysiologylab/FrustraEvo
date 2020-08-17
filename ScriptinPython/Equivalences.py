import sys

pos=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Posiciones",'r')

for linea in pos.readlines():
	if linea[0] == ">":
		salida=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Equivalences/SalidaSRes"+linea[1]+linea[2]+linea[3]+linea[4]+linea[5]+linea[6]+".txt",'w')
		sres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+linea[1]+linea[2]+linea[3]+linea[4]+linea[5]+linea[6]+".pdb.done/FrustrationData/"+linea[1]+linea[2]+linea[3]+linea[4]+linea[5]+linea[6]+".pdb_singleresidue",'r')
		sres.readline()
	else:
		ter=0
		slinea=linea.split(" ")
		tam=len(slinea)
		while(ter<tam-1):
			ter = ter + 1
			if slinea[ter-1] == "G" or slinea[ter-1] == "Z":
				salida.write(str(ter)+"\tN/A\tN/A\tN/A\tN/A\n")
			else:
				sresline=sres.readline()
				sresline = sresline[:-1]
				splitres=sresline.split(" ")
				if splitres[0] != slinea[ter-1]:
					while():
						sres.readline()
						sresline = sresline[:-1]
						splitres=sres.split(" ")
						if splitres[0] != slinea[ter-1] or len(sresline)<1: 
							break
				if splitres[7] > 0.55:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"MIN\n")
				elif splitres[7] < -1:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"MAX\n")
				else:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"NEU\n")
salida.close()
sres.close()

pos.close()
