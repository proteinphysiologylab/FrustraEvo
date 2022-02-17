import sys
import os.path as path
import os

pos=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/ListaPDB.txt",'r')
pathPDB="/home/maria/Documentos/AutoLogo/Pdbs_globins/"
#pathPDB="/media/maria/94e9d876-5628-45e9-96a9-f69bef9cdecc/home/maria/Desktop/PDBs/"
sp=""
frustdir=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"
for PDB in pos.readlines():
	PDB=PDB[:-1]
	sp=PDB.split("_")
	frustra=open(frustdir+"FrustraR.R",'w')
	if sp[0] != ">":	
		tam=len(sp)
		#sp[0]=sp[0].lower()
		if path.exists(pathPDB+sp[0]+".pdb"):
        		cp="cp "+pathPDB+sp[0]+".pdb "+frustdir+sp[0]+".pdb"
			os.system(cp)
     	   	else:
         		wget="cd "+frustdir+";wget \'http://www.rcsb.org/pdb/files/"+sp[0]+".pdb\' -O "+frustdir+sp[0]+".pdb"
         		os.system(wget)
		if tam == 1:
			frustra.write('library(frustratometeR)\nOrderList <- c("'+sp[0]+'.pdb")\nPdb_sr <- dir_frustration(PdbsDir = \"'+frustdir+'\", OrderList = OrderList, Mode = "singleresidue", ResultsDir = \"'+frustdir+'\")')
			PDBid=sp[0]
		else:
			frustra.write('library(frustratometeR)\nOrderList <- c("'+sp[0]+'.pdb")\nPdb_sr <- dir_frustration(PdbsDir = \"'+frustdir+'\", Chain = \"'+sp[1]+'\",OrderList = OrderList, Mode = "singleresidue", ResultsDir = \"'+frustdir+'\")')
			PDBid=sp[0]+'_'+sp[1]
		frustra.close()
		if path.exists(frustdir+sp[0]+".pdb"):
			rsc='cd '+frustdir+';Rscript FrustraR.R'
			os.system(rsc)
		cp="cp "+frustdir+PDBid+".done/FrustrationData/"+PDBid+".pdb_singleresidue "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+PDBid+".done/FrustrationData/"+PDBid+".pdb_msingleresidue"
		os.system(cp)
		rm='cd '+frustdir+';rm Modes.log'
		os.system(rm)
		cp="cp "+frustdir+PDBid+".done/FrustrationData/"+PDBid+".pdb "+sys.argv[1]+"/OutPut"+sys.argv[2]+"/PDB/"+PDBid+".pdb"
		os.system(cp)
pos.close()

#----ChangesInAligment---

fixa="python3 "+sys.argv[1]+"/ScriptinPython/FixAlign.py "+sys.argv[1]+" "+sys.argv[2]
os.system(fixa)
verfrust="python3 "+sys.argv[1]+"/ScriptinPython/FrustraVerif.py "+sys.argv[1]+" "+sys.argv[2]
os.system(verfrust)

#
delete="python3 "+sys.argv[1]+"/ScriptinPython/DeleteGaps.py  "+sys.argv[1]+" "+sys.argv[2]
os.system(delete)

#---------HHMSearch----

hmmbuild = "hmmbuild --amino "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/familias.hmm"+" "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignClean.fasta"

os.system(hmmbuild)
hmmsearch = "hmmsearch --domtblout "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/salida"+" "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/familias.hmm "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignSearch.fasta"
os.system(hmmsearch)

hmmsal = open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/salida",'r')

#-- changes in alignment

ini=0
pdbidsal=""

for line in hmmsal.readlines():
	if line[0] != "#":
		spline= line.split(" ")
		pdbidsal=spline[0]
		break
hmmsal.close()
seqal=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign.fasta",'r')
seq=""
#pdbidsal="1hbh"

while True:
	seqali=seqal.readline()
	if not line:
		break
	seqali=seqali[:-1]
	t=len(seqali)
	pdb = seqali[1:]
	if pdb == pdbidsal:
		seqali=seqal.readline()
		seq=seqali
		break
	seqali=seqal.readline()

seqal.close()

align=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign.fasta",'r')
sal=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign2.fasta",'w')

sal.write(">"+pdbidsal+"\n"+seq)

for line in align.readlines():
	line=line[:-1]
	if line[0] == ">":
		line = line[:1]
		if line  == seq:
			sal.write(line+"\n")
			line=align.readline()
			sal.write(line)

align.close()

align=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign.fasta",'r')

for line in align.readlines():
	sal.write(line)

align.close()
sal.close()

#--- Missing Residuos 

if sys.argv[3] == "Y":
	mc="python2 "+sys.argv[1]+"/ScriptinPython/MissingComplete.py "+sys.argv[1]+"/ "+sys.argv[2]
	os.system(mc)
	ca="python3 "+sys.argv[1]+"/ScriptinPython/ChangeAlign.py "+sys.argv[1]+" "+sys.argv[2]
	os.system(ca)
else:
	cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign2.fasta "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign3.fasta "
	os.system(cp)
	cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign2.fasta "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign4.fasta "
	os.system(cp)
	
