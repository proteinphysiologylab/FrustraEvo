import sys

pos=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Posiciones",'r')

for linea in pos.readlines():
	linea=linea[:-1]
	if linea[0] == '>':
		line=linea[1:]
		spline=line.split('_')
		pdbch=spline[0]
		print(spline[0])
		if len(spline)>=2:
			ch=spline[1]
			pdbch=spline[0]+'_'+spline[1]
		salida=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Equivalences/SalidaSRes"+pdbch+".txt",'w')
		sres=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+pdbch+".pdb.done/FrustrationData/"+spline[0]+".pdb_singleresidue",'r')
		sresline=sres.readline()
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
				if len(splitres) < 7:
					break
				if splitres[0] != slinea[ter-1]:
					while True:
						sresline=sres.readline()
						sresline = sresline[:-1]
						splitres=sresline.split(" ")
						if splitres[0] == slinea[ter-1] or len(sresline)<1:
							print (slinea[ter-1])
							break
				if float(splitres[7]) > 0.55:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"MIN\n")
				elif float(splitres[7]) < -1:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"MAX\n")
				else:
					salida.write(str(ter)+"\t"+str(splitres[0])+"\t"+str(splitres[3])+"\t"+str(splitres[7],)+"\t"+"NEU\n")
salida.close()
sres.close()

pos.close()
