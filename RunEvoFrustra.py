import sys
import os

fastafile=sys.argv[2]
jobID=sys.argv[1]
jobsDir=os.getcwd()

rm='rm -r '+jobsDir+'/OutPutFiles'+jobID
os.system(rm)
rm='rm -r '+jobsDir+'/OutPut'+jobID
os.system(rm)

missa = input('Considers Missing Residues (Y/N): ')
missa=missa.upper()
#missa='N'
d=1

while d:
	if missa == "Y" or missa == "N":
		d=0
	else:
		missa = input('Considers Missing Residues (Y/N): ')
		missa=missa.upper()

#--- Make directories----
pathm='mkdir '+jobsDir+'/OutPut'+jobID
os.system(pathm)
pathm='mkdir '+jobsDir+'/OutPut'+jobID+'/PDB'
os.system(pathm)

pathm='mkdir '+jobsDir+'/OutPutFiles'+jobID
os.system(pathm)
pathm='mkdir '+jobsDir+'/OutPutFiles'+jobID+'/Frustration'
os.system(pathm)
pathm='mkdir '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences'
os.system(pathm)

#--- Copy files ---
cp='cp '+jobsDir+'/ScriptinPython/Generator.R '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/Generator.R'
os.system(cp)
cp='cp '+jobsDir+'/ScriptinPython/SeqLogo.R '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/SeqLogo.R'
os.system(cp)
cp='cp '+jobsDir+'/ScriptinPython/Logo.R '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/Logo.R'
os.system(cp)

c=0
n=0

#--- Changing Alignment --- 

salidafasta=open(jobsDir+'/OutPutFiles'+jobID+'/Sequences.fasta','w') # Sequences without gaps 
sfasta = open(jobsDir+'/'+fastafile,'r') # Input Alignment 
alin = open(jobsDir+'/OutPutFiles'+jobID+'/Alignment.fasta','w') # Alignment without \n
co=0

for linef in sfasta.readlines():
	if linef[0] == '>':
		alin.write('>Seq'+str(c)+'\n')
	else:
		co=0
		lon=len(linef)
		while co < lon:
			if linef[co] != '-':
				salidafasta.write(linef[co])
			alin.write(linef[co])
				
			co+=1
salidafasta.close()
sfasta.close()
alin.close()


tamanio=len(sys.argv)


#--- Run Frustration --- 

frustra='python2 '+jobsDir+'/ScriptinPython/FrustraPDB.py '+jobsDir+' '+jobID+' '+missa
os.system(frustra)

#----Equivalences---

final='python3 '+jobsDir+'/ScriptinPython/FinalAlign.py '+jobsDir+' '+jobID
os.system(final)
equivalences='python3 '+jobsDir+'/ScriptinPython/Equivalences.py '+jobsDir+' '+jobID 
os.system(equivalences)

#-- Call all .R scripts---

fastam='python3 '+jobsDir+'/ScriptinPython/FastaMod.py '+jobsDir+' '+jobID
os.system(fastam)
cd='cd '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences;cat *SalidaSRes* > AllSalidaSResB.txt'
os.system(cd)
logoc='python3 '+jobsDir+'/ScriptinPython/LogoCheck.py '+jobsDir+' '+jobID
os.system(logoc)
logo = 'cd '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences;Rscript Logo.R'
os.system(logo)
gene = 'cd '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences;Rscript Generator.R'
os.system(gene)
#system("cd $jobsDir/OutPutFiles$jobID/Equivalences Rscript SeqLogo.R")

VScript='python3 '+jobsDir+'/ScriptinPython/VScript.py '+jobsDir+' '+jobID #*
os.system(VScript)

cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/HistogramFrustration.svg '+jobsDir+'/OutPut'+jobID+'/HistogramFrustration.svg'
os.system(cp)

cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/SalidaAlign.fasta '+jobsDir+'/OutPut'+jobID+'/SalidaAlign.fasta'
os.system(cp)
cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/EvoFrustra-log.txt '+jobsDir+'/OutPut'+jobID+'/EvoFrustra-log.fasta'
os.system(cp)
cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/ListaPDBC.txt '+jobsDir+'/OutPut'+jobID+'/ListaPDB.txt'
os.system(cp)
cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/CharactPosDataN '+jobsDir+'/OutPut'+jobID+'/IC'
os.system(cp)
cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/Equivalences/IC.csv '+jobsDir+'/OutPut'+jobID+'/IC.csv'
os.system(cp)
cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/familias.hmm '+jobsDir+'/OutPut'+jobID+'/familias.hmm'
os.system(cp)
#cp='cp '+jobsDir+'/OutPutFiles'+jobID+'/seqlogo.png '+jobsDir+'/OutPut'+jobID+'/seqlogo.png
#os.system(cp)

tar='cd '+jobsDir+';tar -zcvf OutPut'+jobID+'.tar.gz OutPut'+jobID
os.system(tar)

#sleep(60)

#system("cd $jobsDir rm -r OutPutFiles$jobID")
#system("cd $jobsDir rm -r OutPut$jobID")
