import sys
import os.path as path
import os

pos=open(sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/ListaPDB.txt",'r')
pathFrustra="/home/maria/bin/frustratometer2-master/"
pathPDB="/home/maria/Descargas/PED1AAA-pdb"
sp=""

for PDB in pos.readlines():
	PDB=PDB[:-1]
	sp=PDB.split("_")
	if sp[0] != ">":	
		tam=len(sp)
		if path.exists("pathPDB/"+sp[0]+".pdb"):
        		cp="cp "+pathPDB+"/"+sp[0]+".pdb "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+sp[0]+".pdb"
			os.system(cp)
     	   	else:
         		wget="wget \'http://www.rcsb.org/pdb/files/"+sp[0]+".pdb\' -O "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+sp[0]+".pdb"
         		os.system(wget)
         		cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+sp[0]+".pdb "+pathPDB+"/"+sp[0]+".pdb"
         		os.system(cp)
		if tam==1:
			cp="cp "+pathPDB+"/"+PDB+".pdb "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+PDB+".pdb"
			os.system(cp)
			cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/"+sp[0]+".pdb "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+".pdb"
			os.system(cp)
			cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+".pdb "+pathFrustra+"/"+sp[0]+".pdb"
			os.system(cp)
			runfrustra="cd "+pathFrustra+"; perl RunFrustratometer.pl "+sp[0]+".pdb singleresidue"
			os.system(runfrustra)
			mvfrustra="cp "+pathFrustra+sp[0]+".pdb.done "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+".pdb.done"
			os.system(mvfrustra)
			mvpdb="cp "+pathFrustra+sp[0]+".pdb "+sys.argv[1]+"/OutPut"+sys.argv[2]+"/PDB/"+sp[0]+".pdb"
			os.system(mvpdb)
			rmall="rm -r "+pathFrustra+sp[0]+".pdb.done"
			os.system(rmall)
			rmall="rm "+pathFrustra+sp[0]+".pdb"
			os.system(rmall)
		else:
			if tam==2:
				sperl="python "+sys.argv[1]+"/ScriptinPython/ChainSeparate.py "+sp[0]+" "+sp[1]+" "+ sys.argv[1]+" "+sys.argv[2]
				os.system(sperl)
			else:
				sperl="python "+sys.argv[1]+"/ScriptinPython/ChainSeparate.py "+sp[0]+" "+sp[1]+" "+ sys.argv[1]+" "+sys.argv[2]+" "+sp[2]
				os.system(sperl)
			cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+".pdb "+pathFrustra+"/"+sp[0]+"_"+sp[1]+".pdb"
			os.system(cp)
			cp="cd "+pathFrustra+";perl RunFrustratometer.pl"+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+".pdb "+pathFrustra+sp[0]+"_"+sp[1]+".pdb"
			os.system(cp)
			runfrustra="cd "+pathFrustra+"; perl RunFrustratometer.pl "+sp[0]+"_"+sp[1]+".pdb singleresidue"
			os.system(runfrustra)	
			mvfrustra="cp -r "+pathFrustra+sp[0]+"_"+sp[1]+".pdb.done "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+"_"+sp[1]+".pdb.done"
			os.system(mvfrustra)
			mvpdb="cp "+pathFrustra+sp[0]+"_"+sp[1]+".pdb "+sys.argv[1]+"/OutPut"+sys.argv[2]+"/PDB/"+sp[0]+".pdb"
			os.system(mvpdb)
			rmall="rm -r "+pathFrustra+sp[0]+"_"+sp[1]+".pdb.done"
			os.system(rmall)
			rmall="rm -r "+pathFrustra+sp[0]+"_"+sp[1]+".pdb"
			os.system(rmall)
			os.system(cp)
			cp="cp "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+"_"+sp[1]+".pdb.done/FrustrationData/"+sp[0]+"_"+sp[1]+".pdb_singleresidue "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/Frustration/"+sp[0]+"_"+sp[1]+".pdb.done/FrustrationData/"+sp[0]+".pdb_singleresidue"
			os.system(cp)
	
pos.close()

#----ChangesInAligment---

fixa="python "+sys.argv[1]+"/ScriptinPython/FixAlign.py "+sys.argv[1]+" "+sys.argv[2]
os.system(fixa)
verfrust="python "+sys.argv[1]+"/ScriptinPython/FrustraVerif.py "+sys.argv[1]+" "+sys.argv[2]
os.system(verfrust)

#
delete="python "+sys.argv[1]+"/ScriptinPython/DeleteGaps.py  "+sys.argv[1]+" "+sys.argv[2]
os.system(delete)

#---------HHMSearch----

hmmbuild = "hmmbuild --amino "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/familias.hmm"+" "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignClean.fasta"

os.system(hmmbuild)
hmmsearch = "hmmsearch --domtblout "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/salida"+" "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/familias.hmm "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/AlignClean.fasta"
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
#my $pdbidsal="1ovm_A_1_186"

for seqali in seqal.readlines():
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

sal.write(">"+pdbidsal+"\n"+seq+"\n")

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
	mc="python "+sys.argv[1]+"/ScriptinPython/MissingComplete.py "+sys.argv[1]+" "+sys.argv[2]
	os.system(mc)
	ca="python "+sys.argv[1]+"/ScriptinPython/ChangeAlign.py "+sys.argv[1]+" "+sys.argv[2]
	os.system(ca)
else:
	cp=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign2.fasta "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign3.fasta "
	os.system(cp)
	cp=sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign2.fasta "+sys.argv[1]+"/OutPutFiles"+sys.argv[2]+"/SeqAlign4.fasta "
	os.system(cp)
	
