#!/usr/bin/python3
import sys
import os.path as path
import os
import math
import numpy as np

##para llenar-----------

#lista de pdbs para analizar
listapdbs=sys.argv[2]+'/PDB_ListChk.txt'
n=int(sys.argv[1])
#archivo de equivalencia para cada pdb (entre numeracion del pdb y posicion del residuo en alineamiento)
equivalenciares='Equival_'

#archivo donde se guardaran los calculos.
archivosalida= 'IC_Conf'

#path donde se encuentran folders con chains frustradas. Las mias eran solo chain A.
frustpath=sys.argv[2]+'/Frustration/'
#path donde se encuentran las equivalencias entre residuos
equivalenciapath=sys.argv[2]+'/Equivalences/'

equival_ref=open(equivalenciapath+equivalenciares+sys.argv[3]+'.txt','r')
lequival=equival_ref.readlines()
equival_ref.close()

matrices=[]

pdbfile=open(equivalenciapath+equivalenciares+sys.argv[3]+'.txt','r')
vect=[]

for line in pdbfile.readlines():
	sp_res=line.split('\t')
	vect.append(sp_res[1])

pdbfile.close()

frstfile=open(frustpath+sys.argv[3]+".done/FrustrationData/"+sys.argv[3]+".pdb_configurational", "r")
pdbfile=open(equivalenciapath+equivalenciares+sys.argv[3].rstrip('\n')+".txt","r")
dic_equival = {}
t_matriz=0

for line in pdbfile.readlines():
	sp_res=line.split('\t')
	dic_equival[sp_res[1]] = sp_res[0]
	if t_matriz < int(sp_res[0]):
		t_matriz = int(sp_res[0])
	if t_matriz < int(sp_res[1]):
		t_matriz = int(sp_res[1])
	
pdbfile.close()

matriz = np.zeros((len(dic_equival)+2,len(dic_equival)+2))
for mi in range(len(matriz)):
	for mj in range(len(matriz)):
		matriz[mi][mj] = -100

for line in frstfile.readlines():
	sp=line.rstrip('\n').split()
	if sp[0] in dic_equival and sp[1] in dic_equival and 'Res' not in line:
		matriz[int(dic_equival[sp[0]])][int(dic_equival[sp[1]])] = float(sp[11])
	
frstfile.close()
matrices.append(matriz)

## funciones ------------

#busca numero de residuo en pdb para cada posicion del alineamiento, ya que el archivo de frustra tiene numeracion pdb.
def read_X(doc,numres,chain): 
	for a in doc:
		lin=a.rstrip('\n')
		lin=a.split()
		#print(ambler,numres,lin[4])
		if int(lin[0]) == int(numres):
			return lin[1]

#busca linea para contacto en frustra file, y suma 1 a la categoria que corresponda (MAX, MIN o NEU)
def cont_X(doc,i,j):
	cut_MIN=0.78
	cut_MAX=-1.0
	numeroconts.append(1)
	frst.append(float(doc))
	if float(doc)>= float(cut_MIN):
		MIN.append(1)
	elif float(doc)<= float(cut_MAX):
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

salida=open(archivosalida, 'w')
salida.write('Res\tRes\tNumeroConts\tFreqConts\tpNEU\tpMIN\tpMAX\tHNEU\tHMIN\tHMAX\tHtotal\tICNEU\tICMIN\tICMAX\tICtotal\tEstadoConservado\n')

for h in range(0,len(listpdb)):
	mat_ref=matrices[0]
	if listpdb[h] != sys.argv[3]:
		pdbfile=open(equivalenciapath+equivalenciares+listpdb[h].rstrip('\n')+".txt","r")
		dic_equival = {}
		t_matriz=0
		for line in pdbfile.readlines():
			if 'N/A' not in line:
				sp_res=line.split('\t')
				dic_equival[sp_res[1]] = sp_res[0]
				if t_matriz < int(sp_res[0]):
					t_matriz = int(sp_res[0])
				if t_matriz < int(sp_res[1]):
					t_matriz = int(sp_res[1])
		pdbfile.close()
	
		frstfile=open(frustpath+listpdb[h].rstrip('\n')+".done/FrustrationData/"+listpdb[h].rstrip('\n')+".pdb_configurational", "r")
		matriz = np.zeros((n+2,n+2))
		for mi in range(0,n+2):
			for mj in range(0,n+2):
				matriz[mi][mj] = -100
	
		for line in frstfile.readlines():
			sp=line.rstrip('\n').split()
			if sp[0] in dic_equival and sp[1] in dic_equival:
				if mat_ref[int(dic_equival[sp[0]])][int(dic_equival[sp[1]])] != -100:
					matriz[int(int(dic_equival[sp[0]]))][int(int(dic_equival[sp[1]]))] = float(sp[11])
		
		frstfile.close()
		matrices.append(matriz)
	
for ri in range(0,len(lequival)):
	lresi=lequival[ri].rstrip('\n')
	resi=lresi.split()
	for rj in range(ri+1,len(lequival)):
		lresj=lequival[rj].rstrip('\n')
		resj=lresj.split()
		numeroconts=[]
		frst=[]
		NEU=[]
		MIN=[]
		MAX=[]		

		for j in range(0,len(listpdb)): #abre cada renglon del archivo, lee nombre pdb, abre pdb y lo lee
			m=matrices[j]
			if m[int(resi[0])][int(resj[0])] != float(-100):
				cont_X(m[int(resi[0])][int(resj[0])],resi[0],resj[0])
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



			salida.write(str(resi[0])+"\t"+str(resj[0])+'\t'+str(conts)+"\t"+str(freqconts)+"\t"+str(pNEU)+'\t'+str(pMIN)+'\t'+str(pMAX)+'\t'+str(HNEU)+'\t'+str(HMIN)+'\t'+str(HMAX)+'\t'+str(Htotal)+'\t'+str(IC_NEU)+'\t'+str(IC_MIN)+'\t'+str(IC_MAX)+'\t'+str(IC_total)+'\t'+str(conservedstate)+"\n")

salida.close()
