#!/usr/bin/python3
import os
import re
import sys
from Bio import SeqIO
import numpy as np


def copyfiles(JodID,path_to_r,path_to_Pdbs):

	"""     This function is for the creation of folders and copying of necessary files
		Parameters:
			- JodID: the job name
			- path_to_r: path to the R files
	"""
	path_direc='FrustraEvo_'+JodID
	if os.path.isdir(path_direc):
		os.system('rm -r '+path_direc)
	os.system('mkdir '+path_direc)
	os.system('mkdir '+path_direc+'/Equivalences')
	os.system('mkdir '+path_direc+'/OutPutFiles')
	os.system('mkdir '+path_direc+'/Frustration')
	os.system('cp '+path_to_r+'/*.R* '+path_direc+'/Equivalences/')
	
def pdb_list(fasta_file):
	"""     This function is for the creation of folders and copying of necessary files
		Parameters:
			- fasta_file: MSA of input
			- pathPDB: path to the folder with the pdbs ()
	"""
	msa=open(fasta_file)
	sp=fasta_file.split('.')
	out=open(sp[0]+'.list','w')
	for lmsa in msa.readlines():
		if lmsa[0] == '>':
			out.write(lmsa[1:])
	out.close()
	msa.close()
	return(sp[0]+'.list')

def changes(JodID,MSA_File):
	"""     This function make changes in the MSA and make a List
		Parameters:
			- JodID: the job name
			- MSA_File: the MSA
			
	"""
	path_direc='FrustraEvo_'+JodID
	out=open(path_direc+'/AlignSE.fasta','w')
	out_list=open(path_direc+'/ListaPDB.txt','w')
	pathAlign=MSA_File
	for seq_record in SeqIO.parse(pathAlign, "fasta"):
		seqid=seq_record.id
		seq=seq_record.seq
		out.write('>'+seqid+'\n'+str(seq)+'\n')
		
	out.close()
	out_list.close()

def FrustraPDB(list_pdbs,JodID,pathPDB='None'):
	"""     This function is for the frustration calculation
		Parameters:
			- list_pdbs: the list with the Pdbs 
			- JodID: the job name
			- pathPDB: path to the folder with the pdbs ()
			
	If you do not have the pdbs structures, the pipeline will download them from the Protein Data Bank database
	"""
	path_direc='FrustraEvo_'+JodID
	pdbs=open(list_pdbs,'r')
	frustdir=path_direc+'/Frustration/'
	for line in pdbs.readlines():
		line=line.rstrip('\n')
		if pathPDB == 'None':
			os.system('cd '+frustdir+'/;wget \'http://www.rcsb.org/pdb/files/'+line+'.pdb\' -O '+JodID+'/Frustration/'+line+'.pdb')
		else:
			#line=line.lower()
			os.system('cp '+pathPDB+'/'+line+'.pdb '+frustdir+'/'+line+'.pdb')
	pdbs.close()
	directory=os.getcwd()+'/'
	frustra=open(frustdir+"FrustraR.R",'w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \"'+directory+frustdir+'\"\nResultsDir <- \"'+directory+frustdir+'\"\ndir_frustration(PdbsDir = PdbsDir, Mode = \"singleresidue\", ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')

  
def checks(JodID):
	"""     This function checks the frustration calculations
	Parameters:
		- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	MSA=open(path_direc+'/AlignSE.fasta','r')
	out=open(path_direc+'/SeqAlign.fasta','w')
	out_list=open(path_direc+'/ListaPDBC.txt','w')	
	out_error=open(path_direc+'/EvoFrustra-log','w')
	error=0
	for line in MSA.readlines():
		line=line.rstrip('\n')
		if line[0] == '>':
			error=0
			line=line[1:]
			path_frst=path_direc+'/Frustration/'+line+'.done/FrustrationData/'+line+'.pdb_singleresidue'
			if os.path.isfile(path_frst):
				frst=open(path_frst,'r')
				cl=0
				for lfrst in frst.readlines():
					cl+=1
					if cl>3:
						error=1
						break
				frst.close()	
			if error == 0:
				out_error.write('>'+line+'\n')
			else:
				out.write('>'+line+'\n')
				out_list.write(line+'\n')		
		elif error==1:
			out.write(line+'\n')
		
	MSA.close()
	out.close()
	out_list.close()
	out_error.close()
	

def DeleteGaps(JodID):
	"""     This function find the reference protein using HHMSearch
		Parameters:
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	AClean=open(path_direc+'/AlignClean.fasta','w')
	Asearch=open(path_direc+'/AlignSearch.fasta','w')
	pathAlign=path_direc+'/SeqAlign.fasta'
	vector = list()
	c=0
	for seq_record in SeqIO.parse(pathAlign, "fasta"):
		c = c + 1
		vector.append(seq_record.seq)
	porc=c*0.6 
	cgaps=0
	aux=""
	vectorgaps = np.zeros(len(vector[0]))
	for i in range(0,len(vector[0])):
		cgaps=0
		for j in range(0,c):
			aux=vector[j]
			if aux[i] == '-':
				cgaps+=1
		if cgaps < porc:
			vectorgaps[i] = 0
		else:
			vectorgaps[i] = 1
	
	for seq_record in SeqIO.parse(pathAlign, "fasta"):
		seqid=seq_record.id
		seq=seq_record.seq
		pr=0
		AClean.write(">"+seqid+"\n")
		Asearch.write(">"+seqid+"\n")
		for i in range(0,len(vectorgaps)):
			if vectorgaps[i]==0:
				if pr==60:
					AClean.write("\n")
					pr=0
				AClean.write(seq[i])
				if seq[i] != '-':
					Asearch.write(seq[i])
				pr+=1
		AClean.write("\n")
		Asearch.write("\n")
	AClean.close()
	Asearch.close()
	os.system('hmmbuild --amino '+path_direc+'/familias.hmm'+' '+path_direc+'/AlignClean.fasta')
	os.system('hmmsearch --domtblout '+path_direc+'/salida'+' '+path_direc+'/familias.hmm '+path_direc+'/AlignSearch.fasta')

	hmmsal = open(path_direc+'/salida','r')
	ini=0
	pdbidsal=""
	for line in hmmsal.readlines():
		if line[0] != "#":
			spline = line.split(" ")
			pdbidsal = spline[0]
			break
	hmmsal.close()
	return pdbidsal


def prepare_file(JodID,prot_ref):
	"""     This function prepare files for the calculations of the Equivalences
	Parameters:
		- JodID: the job name
		- prot_ref: the reference protein
	"""
	path_direc='FrustraEvo_'+JodID
	MSA=open(path_direc+'/SeqAlign.fasta','r')
	out=open(path_direc+'/SeqAlign_Ref.fasta','w')
	f=0
	for line in MSA.readlines():
		line=line.rstrip('\n')
		if line[0] == '>':
			seq=line[1:]
			if seq == prot_ref:
				f=1
				out.write(line+'\n')
		elif f == 1:
			out.write(line+'\n')
			break
	MSA.close()
	out.close()
	os.system('cat '+path_direc+'/SeqAlign.fasta >> '+path_direc+'/SeqAlign_Ref.fasta')


def MissingRes(JodID):
	"""     This function find Missing Residues
		Parameters:
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	hash = {'CYS' : 'C', 'ASP' : 'D', 'SER' : 'S','GLN' : 'Q','LYS' : 'K','ILE' : 'I','PRO' : 'P','THR' : 'T','PHE' : 'F','ASN' : 'N','GLY' : 'G','HIS' : 'H','LEU' : 'L','ARG' : 'R','TRP' : 'W','ALA' : 'A','VAL' : 'V','GLU' : 'E','TYR' : 'Y','MET' : 'M','MSE' : 'B'}
	MSA=open(path_direc+'/SeqAlign_Ref.fasta','r')
	for line in MSA.readlines():
		line=line.rstrip('\n')
		if line[0] == '>':
			pdbch=line[1:]
			spline=pdbch.split('-')
			count=0
			bus=open(JodID+'/Frustration/'+spline[0]+'.pdb','r')
			for busca in bus.readlines():
				spbusca=busca.split(' ')
				if len(spline) >= 2:
					if spbusca[0] == 'REMARK' and spbusca[1]=='465' and spbusca[7] == spline[1]:
						chain.append(ch)
						posi=busca[21]+busca[22]+busca[23]+busca[24]+busca[25]
						com.append(posi.lstrip())
						iden.append(hash[spbusca[6]])
						count=count+1
				else:
					if spbusca[0] == 'REMARK' and spbusca[1]=='465':
						com.append(spbusca[10])
						chain.append(spbusca[7])
						count=count+1
			bus.close()
			if count!=0:
				os.system('cp '+JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue '+JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue')
				sres=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue','r')
				sressal=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','w')
				cnt=0
				spl=''
				cnt2=0
				Sres=sres.readline()
				sressal.write('#Res ChainRes DensityRes AA NativeEnergy DecoyEnergy SDEnergy FrstIndex\n')
				for Sres in sres.readlines():
					spl= Sres.split()
					if spl[0] == '#Res':
						sressal.write(Sres)
					else:
						if spl[0] > com[cnt]:
							while (int(spl[0]) - 1 >= int(com[cnt]) and count>cnt2):
								sressal.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n')
								cnt=cnt+1
								cnt2=cnt2+1
								if cnt == count:
									cnt=cnt-1
									break
							
							sressal.write(Sres)
						else:
							sressal.write(Sres)
				if spl[0]<com[cnt]:
					while cnt<=count-1:
						sressal.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n')
						cnt=cnt+1
				sres.close()
				sressal.close()
	align.close()

def ChangeAlign():
	"""     This function make changes in the alignment in the case if you hace missing residues
	Parameters: 
		- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	os.system('cp '+path_direc+'/SeqAlign_Ref.fasta '+path_direc+'/SeqAlign_Ref_aux.fasta')
	pathAlign=path_direc+'/SeqAlign_Ref_aux.fasta'

	archisal3=open(path_direc+'/SeqAlign_Ref_aux_1.fasta','w')
	archisal4=open(path_direc+'/SeqAlign_Ref.fasta','w')

	for seq_record in SeqIO.parse(pathAlign, "fasta"):
		seqid=seq_record.id
		seq=seq_record.seq
		pdbid=seqid
		archisal3.write(">"+seqid+"\n")	
		archisal4.write(">"+seqid+"\n")
		d=len(sseqid)
		archisres=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r')
		c=0
		if(d>2):
			c=int(sseqid[d-2])
		archisres.readline()
		l = archisres.readline()
		sl = l.split(" ")
		i=int(sl[0])
		while int(i)<int(c):
			l = archisres.readline()
			i = i + 1
		w=0
		while True:
			srespos =  l.split(" ") 
			if seq[w] == '-' or seq[w] == 'X' or seq[w] == 'B':
				archisal3.write("-")
				archisal4.write("-")
			elif len(srespos) >= 4 and srespos[4] == "Missing":
					archisal3.write("-")
					archisal4.write("Z")
					l = archisres.readline()
			else:
				archisal3.write(seq[w])
				archisal4.write(seq[w])
				l = archisres.readline()				
			w = w + 1		

			if w == len(seq):	
				break
		archisal4.write("\n")
		archisal3.write("\n")
	archisal3.close()
	archisal4.close()
	archisres.close()
	os.system('rm '+path_direc+'/SeqAlign_Ref_aux_1.fasta')
	os.system('rm '+path_direc+'/SeqAlign_Ref_aux.fasta')
	

def FinalAlign(JodID):
	"""     This function create the final files
	Parameters: 
		- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	pathAlign=path_direc+'/SeqAlign_Ref.fasta'
	salida1=open(path_direc+'/SalidaAlign.fasta','w')
	salida2=open(path_direc+'/SalidaAlignSE.fasta','w')
	salida3=open(path_direc+'/Posiciones','w')

	vector = list()
	c=0	
	m=0
	control=0 
	cgaps=0;
	tam=0
	tam=0
	r=0

	for seq_record in SeqIO.parse(pathAlign, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq
		pdb=seqid
		sres=''
		longline=0
		if c!=0:
			salida1.write('>'+pdb+'\n')
			salida2.write('>'+pdb+'\n')
			salida3.write('>'+pdb+'\n')
			sres=open(path_direc+'/Frustration/'+pdb+'.done/FrustrationData/'+pdb+'.pdb_singleresidue','r')
			lsres=sres.readline()
		if c==0:
			c=1
			tam=len(seq)
			i=0
			while(i<tam):
				if seq[i] == '-' or seq[i] == 'Z':
					vector.append(0)
				else:
					vector.append(1)
				i += 1
		else:
			splitres=''
			tam=len(seq)
			lsres = sres.readline()
			if r==1 and len(pdbid)>2:
				splitres=lsres.split(' ')
				r= int(splitres[0])
				a = int(pdbid[2])-1
				while (r<a):
					lsres=sres.readline()
					r = r + 1
			n=0
			q=0
			for j in range (0,tam):
				if vector[j] == 0:
					if seq[j] != '-':
						lsres = sres.readline()
				else:
					if seq[j] == 'Z':
						salida1.write('-')
						longline+=1
						n = n + 1
						q = q + 1
						salida2.write('-')
					else:
						salida1.write(str(seq[j]))
						longline+=1
						n = n + 1
						q = q + 1
						salida2.write(str(seq[j]))
					if n==51:
						salida1.write('\n')
						n=0
					if seq[j] == 'Z':
						salida3.write('Z ')
					elif seq[j] == '-':
						salida3.write('G ')
					else:
						if len(splitres) < 2:
							splitres = lsres.split(' ')
						lsres = sres.readline()
						salida3.write(str(splitres[0])+' ')
						splitres = lsres.split(' ')
						r = r + 1
						
			salida1.write('\n')
			salida2.write('\n')
			salida3.write('\n')
	sres.close()


	longalign=open(path_direc+'/Equivalences/long.txt','w')
	longalign.write(str(longline)+'\n')
	longalign.close()
	os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/long.txt')

	salida1.close()
	salida2.close()
	salida3.close()
	

def Equivalences(JodID):
	"""     This function create the final files
		Parameters: 
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	pos=open(path_direc+'/Posiciones','r')

	for linea in pos.readlines():
		linea=linea[:-1]
		if linea[0] == '>':
			line=linea[1:]
			pdbch=line
			salida=open(path_direc+'/Equivalences/SalidaSRes'+pdbch+'.txt','w')
			sres=open(path_direc+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r')
			sresline=sres.readline()
		else:
			ter=0
			slinea=linea.split(" ")
			splitres=sresline.split(" ")
			tam=len(slinea)
			while(ter<tam - 1):
				ter+=1
				if slinea[ter-1] == "G" or slinea[ter-1] == "Z":
					salida.write(str(ter)+'\tN/A\tN/A\tN/A\tN/A '+splitres[1]+'\n')
					
				else:	
					sresline=sres.readline()
					sresline = sresline.rstrip('\n')
					splitres=sresline.split(" ")
					if sresline != '':
						if len(splitres) < 7 and splitres[4] !='Missing':
							break
						if int(splitres[0]) < int(slinea[ter-1]):
							while True:
								sresline=sres.readline()
								sresline = sresline.rstrip('\n')
								splitres=sresline.split(" ")
								if splitres[0] == slinea[ter-1] or len(sresline)<1:
									break
						if len(splitres) == 6:
							ter-=1
						if splitres[0] == slinea[ter-1] and len(splitres) > 7:
							if float(splitres[7]) > 0.55:
								salida.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MIN '+splitres[1]+'\n')
							elif float(splitres[7]) < -1:
								salida.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MAX '+splitres[1]+'\n')
							else:
								salida.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'NEU '+splitres[1]+'\n')
	salida.close()
	sres.close()

	pos.close()
	
	
def FastaMod(JodID):
	"""     This function create the file for the sequence logo: 
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	fasta=path_direc+'/SalidaAlignSE.fasta'
	archisal=open(path_direc+'/Equivalences/Logo.fasta','w')

	for seq_record in SeqIO.parse(fasta, "fasta"):
		seqid=seq_record.id
		seq=seq_record.seq
		archisal.write(str(seq)+"\n")
	archisal.close()


def LogoCheck(JodID):
	"""     This function checks the logo 
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	os.system('cd '+path_direc+'/Equivalences/;cat *SalidaSRes* > AllSalidaSResB.txt')
	res=open(path_direc+'/Equivalences/AllSalidaSResB.txt','r')
	sal=open(path_direc+'/Equivalences/AllSalidaSRes.txt','w')
	for lres in res:
		lres = lres[:-1]
		splres= lres.split('\t')
		tam=len(splres)
		if tam>4:
			sal.write(str(lres)+'\n')

	res.close()
	sal.close()


def plots_logo(JodID,prot_ref):
	"""     This function makes the plots for the frustration logo: 
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	os.system('cd '+path_direc+'/Equivalences;Rscript Logo.R')
	add_ref(JodID,prot_ref)
	os.system('cd '+path_direc+'/Equivalences;Rscript Generator.R')
	os.system('mv '+path_direc+'/Equivalences/HistogramFrustration.svg'+' '+path_direc+'/OutPutFiles/FrustrationLogo'+JodID+'.svg')
	#os.system('tar -zcvf '+JodID+'.tar.gz '+path_direc)

def add_ref(JodID,prot_ref):
	"""     This function generate the .pml files: 
			- JodID: the job name
			- prot_ref: the reference protein
	"""
	path_direc='FrustraEvo_'+JodID
	res=open(path_direc+'/Equivalences/SalidaSRes'+prot_ref+'.txt','r')
	out=open(path_direc+'/Equivalences/CharactPosData','w')
	ic=open(path_direc+'/Equivalences/CharactPosDataN','r')
	vectorAA=[]
	num=[]
	for line in res.readlines():
		line=line.rstrip('\n')
		sp=line.split()
		vectorAA.append(sp[2])
		num.append(sp[1])
	out.write('Res	AA_Ref	Num_Ref	Prot_Ref	%Min	%Neu	%Max	CantMin	CantNeu	CantMax	ICMin	ICNeu	ICMax	ICTot 	 FrustEstado\n')
	a=0
	for line in ic.readlines():
		line=line.rstrip('\n')
		if a==0:
			a+=1
		else:
			sp=line.split()
			out.write(sp[0]+'\t'+vectorAA[int(sp[0])-1]+'\t'+str(num[int(sp[0])-1])+'\t'+prot_ref+'\t'+sp[1]+'\t'+sp[2]+'\t'+sp[3]+'\t'+sp[4]+'\t'+sp[5]+'\t'+sp[6]+'\t'+sp[7]+'\t'+sp[8]+'\t'+sp[9]+'\t'+sp[10]+'\t'+sp[11]+'\n')

	ic.close()
	out.close()	
	res.close()
	os.system('cp '+path_direc+'/Equivalences/CharactPosData'+' '+path_direc+'/OutPutFiles/IC_SingleRes_'+JodID)

def VScript(JodID):
	"""     This function generate the .pml files: 
			- JodID: the job name
	"""
	path_direc='FrustraEvo_'+JodID
	Lista= open(path_direc+'/ListaPDBC.txt','r')
	os.system('mkdir '+path_direc+'/Visualization')
	os.system('cp '+path_direc+'/Frustration/*.pdb* '+path_direc+'/Visualization/')
	for lista in Lista:
		lista=lista[:-1]
		pdbid=lista
		sal = open(path_direc+'/Visualization/'+pdbid+'.pml','w')
		Equ = open(path_direc+'/Equivalences/SalidaSRes'+pdbid+'.txt','r')
		ECon = open(path_direc+'/Equivalences/CharactPosDataN','r')
		sal.write('load '+pdbid+'.pdb\nhide all\nshow cartoon, all\nbg_color white\ncolor black, all')
		EstCon=ECon.readline()
		for EstCon in ECon.readlines():
			Equi= Equ.readline()
			splitE= Equi.split()
			splitEC= EstCon.split()
			if splitEC[11] == 'MAX' and float(splitEC[10])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor red,resi '+splitE[1]+' and chain '+splitE[5])
					
			if splitEC[11] == 'MIN' and float(splitEC[10])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor green,resi '+splitE[1]+' and chain '+splitE[5])

			if splitEC[11] == 'NEU' and float(splitEC[10])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor gray,resi '+splitE[1]+' and chain '+splitE[5])
		sal.close()
		Equ.close()
		ECon.close()
	Lista.close()
	
def CMaps_Mutational(JodID,path_to_r,prot_ref):
	"""     This function generate the CMaps for mutational index: 
			- JodID: the job name
			- path_to_r: path to the R files
	"""
	directory=os.getcwd()+'/'
	path_direc='FrustraEvo_'+JodID
	frustdir='FrustraEvo_'+JodID+'/Frustration/'
	if os.path.isdir(path_direc+'/CMaps'):
		print('The Path: '+path_direc+'/CMaps exists')
	else:
		os.system('mkdir '+path_direc+'/CMaps')
	os.system('cp '+path_direc+'/ListaPDBC.txt '+path_direc+'/CMaps/ListaPDBC.txt')
	os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
	os.system('cp '+path_to_r+'/*.r* '+path_direc+'/CMaps/')
	os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
	frustra=open(frustdir+"FrustraR.R",'w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \"'+directory+frustdir+'\"\nResultsDir <- \"'+directory+frustdir+'\"\ndir_frustration(PdbsDir = PdbsDir, Mode = \"mutational\", ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
	longalign=open(path_direc+'/Equivalences/long.txt','r')
	long=longalign.readline()
	long=long[:-1]
	longalign.close()
	dir_total=directory+path_direc
	os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Mutational.py '+str(long)+' '+dir_total)
	os.system('cd '+path_direc+'/CMaps;Rscript IC_conts_Mut.r')
	add_ref_Cmaps(JodID,prot_ref,'Mut')
	os.system('cp '+path_direc+'/CMaps/IC_Mut.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JodID+'_Mut.png')

	
def CMaps_Configurational(JodID,path_to_r,prot_ref):
	"""     This function generate the CMaps for Configurational index: 
			- JodID: the job name
			- path_to_r: path to the R files
	"""
	directory=os.getcwd()+'/'
	path_direc='FrustraEvo_'+JodID
	frustdir='FrustraEvo_'+JodID+'/Frustration/'
	if os.path.isdir(path_direc+'/CMaps'):
		print('The Path: '+path_direc+'/CMaps exists')
	else:
		os.system('mkdir '+path_direc+'/CMaps')
		os.system('cp '+path_direc+'/ListaPDBC.txt '+path_direc+'/CMaps/ListaPDBC.txt')
		os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
		os.system('cp '+path_to_r+'/*.r* '+path_direc+'/CMaps/')
		os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
	frustra=open(frustdir+"FrustraR.R",'w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \"'+directory+frustdir+'\"\nResultsDir <- \"'+directory+frustdir+'\"\ndir_frustration(PdbsDir = PdbsDir, Mode = \"configurational\", ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
	longalign=open(path_direc+'/Equivalences/long.txt','r')
	long=longalign.readline()
	long=long[:-1]
	longalign.close()
	dir_total=directory+path_direc
	os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Conf.py '+str(long)+' '+dir_total)
	os.system('cd '+path_direc+'/CMaps;Rscript IC_conts_Conf.r')
	os.system('cp '+path_direc+'/CMaps/IC_Conf.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JodID+'_Conf.png')
	add_ref_Cmaps(JodID,prot_ref,'Conf')
	
def add_ref_Cmaps(JodID,prot_ref,Mode):
	"""     This function generate the .pml files: 
			- JodID: the job name
			- prot_ref: the reference protein
	"""
	path_direc='FrustraEvo_'+JodID
	res=open(path_direc+'/Equivalences/SalidaSRes'+prot_ref+'.txt','r')
	out=open(path_direc+'/CMaps/IC_'+Mode+'_ref','w')
	ic=open(path_direc+'/CMaps/IC_'+Mode,'r')
	vectorAA=[]
	num=[]
	for line in res.readlines():
		line=line.rstrip('\n')
		sp=line.split()
		vectorAA.append(sp[2])
		num.append(sp[1])
	out.write('Res	Res AA1	AA2	NumRes1_Ref	NumRes2_Ref	Prot_Ref	NumeroConts	FreqConts	pNEU	pMIN	pMAX	HNEU	HMIN	HMAX	Htotal	ICNEU	ICMIN	ICMAX	ICtotal	EstadoConservado\n')
	a=0
	for line in ic.readlines():
		line=line.rstrip('\n')
		if a==0:
			a+=1
		else:
			sp=line.split()
			out.write(sp[0]+'\t'+sp[1]+'\t'+vectorAA[int(sp[0])-1]+'\t'+vectorAA[int(sp[1])-1]+'\t'+str(num[int(sp[0])-1])+'\t'+str(num[int(sp[1])-1])+'\t'+prot_ref+'\t'+sp[2]+'\t'+sp[3]+'\t'+sp[4]+'\t'+sp[5]+'\t'+sp[6]+'\t'+sp[7]+'\t'+sp[8]+'\t'+sp[9]+'\t'+sp[10]+'\t'+sp[11]+'\t'+sp[12]+'\t'+sp[13]+'\t'+sp[14]+'\t'+sp[15]+'\n')

	ic.close()
	out.close()	
	res.close()
	os.system('cp '+path_direc+'/CMaps/IC_'+Mode+'_ref'+' '+path_direc+'/OutPutFiles/IC_'+Mode+'_'+JodID)
