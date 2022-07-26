#!/usr/bin/python3
import sys
import os.path as path
import os
import math

##para llenar-----------

#rango de residuos
posicionmin=1
posicionmax=int(sys.argv[1])+1

#lista de pdbs para analizar
listapdbs='PDB_ListChk.txt'

#archivo de equivalencia para cada pdb (entre numeracion del pdb y posicion del residuo en alineamiento)
equivalenciares='Equival_'

#archivo donde se guardaran los calculos.
archivosalida= 'IC_Conf'

#path donde se encuentran folders con chains frustradas. Las mias eran solo chain A.
frustpath=sys.argv[2]+'/Frustration/'
#path donde se encuentran las equivalencias entre residuos
equivalenciapath=sys.argv[2]+'/Equivalences/'


## funciones ------------

#busca numero de residuo en pdb para cada posicion del alineamiento, ya que el archivo de frustra tiene numeracion pdb.
def read_X(doc,numres): 
	for a in doc:
		lin=a.split("\t")
		ambler=str(lin[0]).rstrip('\n')
		#print(ambler,numres)
		if str(ambler) == str(numres):
			return str(lin[1])

#busca linea para contacto en frustra file, y suma 1 a la categoria que corresponda (MAX, MIN o NEU)
def cont_X(doc,i,j):
	cut_MIN=0.78
	cut_MAX=-1.0
	for a in doc:
		lin=a.split(" ")
		#print(lin[0], i, lin[1], j, lin[11])
		if str(lin[0])==str(i) and str(lin[1])==str(j):		
			numeroconts.append(1)
#si dice lin[2], calculo IC de configurational entropy. Si dice lin[3] es de mutational entropy	
			frst.append(float(lin[11]))
			if float(lin[11])>= float(cut_MIN):
				MIN.append(1)
			elif float(lin[11])<= float(cut_MAX):
				MAX.append(1)
			else:
				NEU.append(1)
			return True

#calculos Shannon 

def Hmin(sumaMIN, pMIN):
	if pMIN>0:
		logMIN=math.log(pMIN, 2)
		HMIN=pMIN*logMIN
		HMIN=-(HMIN)
	else:
		HMIN=0
	return HMIN

def Hmax(sumaMAX, pMAX):
	if pMAX>0:
		logMAX=math.log(pMAX, 2)
		HMAX=pMAX*logMAX
		HMAX=-(HMAX)
	else:
		HMAX=0
	return HMAX

def Hneu(sumaNEU, pNEU):
	if pNEU>0:
		logNEU=math.log(pNEU, 2)
		HNEU=pNEU*logNEU
		HNEU=-(HNEU)
	else:
		HNEU=0
	return HNEU
	
def information_content(Htot):
	pmin_esperada=0.4
	logmin=math.log(pmin_esperada, 2)
	pmax_esperada=0.1
	logmax=math.log(pmax_esperada, 2)
	pneu_esperada=0.5
	logneu=math.log(pneu_esperada, 2)
	Hmin_esperada= pmin_esperada * logmin
	Hmax_esperada= pmax_esperada * logmax
	Hneu_esperada= pneu_esperada * logneu
	Htot_esperada= -(Hmin_esperada + Hmax_esperada + Hneu_esperada)
	informationcontent=Htot_esperada- Htot
	return informationcontent



## --------------------------

listapdbs=open(listapdbs,"r")
listpdb=listapdbs.readlines()
listapdbs.close()
total=len(listpdb)

salida=open(archivosalida, 'a')
salida.write('Res\tRes\tNumeroConts\tFreqConts\tpNEU\tpMIN\tpMAX\tHNEU\tHMIN\tHMAX\tHtotal\tICNEU\tICMIN\tICMAX\tICtotal\tEstadoConservado\n')

for resi in range(posicionmin, posicionmax):
	for resj in range(posicionmin, posicionmax):
#pongo contadores
		numeroconts=[]
		frst=[]
		NEU=[]
		MIN=[]
		MAX=[]


		for j in listpdb: #abre cada renglon del archivo, lee nombre pdb, abre pdb y lo lee
			pdb=j.rstrip('\n')			
			pdbfile=open(equivalenciapath+equivalenciares+pdb+".txt","r")
			arch=pdbfile.readlines() #crea lista cn lineas del doc
			pdbfile.close()
			ipdb=read_X(arch,resi)
			jpdb=read_X(arch,resj)
			if ipdb!=None and jpdb!=None: #si ese res ambler esta en la estruct
				frstfile=open(frustpath+pdb+".done/FrustrationData/"+pdb+".pdb_configurational", "r")
				doc=frstfile.readlines()[1:]
				frstfile.close()
				cont_X(doc,ipdb,jpdb)
		conts=sum(numeroconts)
		num_MIN=sum(MIN)
		num_MAX=sum(MAX)
		num_NEU=sum(NEU)

	#si el contacto aparece una sola vez, no tiene sentido calcular h e ic.. 
		if conts>1:			
		
			if max(num_MIN, num_NEU, num_MAX)==num_MIN:
				conservedstate='MIN'
			elif max(num_MIN, num_NEU,num_MAX)==num_NEU:
				conservedstate='NEU'
			elif max(num_MIN, num_NEU,num_MAX)==num_MAX:
				conservedstate='MAX'

			
			pNEU=float(num_NEU)/float(conts)
			pMIN=float(num_MIN)/float(conts)
			pMAX=float(num_MAX)/float(conts)
			
			HMIN=Hmin(num_MIN, pMIN)
			HMAX=Hmax(num_MAX, pMAX)
			HNEU=Hneu(num_NEU, pNEU)
			Htotal=HMIN+HMAX+HNEU

			IC_total=information_content(Htotal)
			IC_MIN=IC_total*pMIN
			IC_MAX=IC_total*pMAX
			IC_NEU=IC_total*pNEU

			freqconts=float(conts)/float(total)



			salida.write(str(resi)+"\t"+str(resj)+'\t'+str(conts)+"\t"+str(freqconts)+"\t"+str(pNEU)+'\t'+str(pMIN)+'\t'+str(pMAX)+'\t'+str(HNEU)+'\t'+str(HMIN)+'\t'+str(HMAX)+'\t'+str(Htotal)+'\t'+str(IC_NEU)+'\t'+str(IC_MIN)+'\t'+str(IC_MAX)+'\t'+str(IC_total)+'\t'+str(conservedstate)+"\n")

salida.close()

