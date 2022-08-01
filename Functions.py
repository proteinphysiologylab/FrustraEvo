#!/usr/bin/python3
import os
import re
import sys
from Bio import SeqIO
import numpy as np


def PreparingFiles(JobId,path_to_r,path_to_Pdbs,fasta):
	copyfiles(JobId,path_to_r,path_to_Pdbs)
	list_pdb=pdb_list(fasta)
	changes(JobId,fasta)
	return list_pdb
	
def copyfiles(JobID,path_to_r,path_to_Pdbs):

	'''     This function is for the creation of folders and copying of necessary files
		Parameters:
			- JobID: the job name
			- path_to_r: path to the R files
	'''
	path_direc='FrustraEvo_'+JobID
	if os.path.isdir(path_direc):
		os.system('rm -r '+path_direc)
	os.system('mkdir '+path_direc)
	os.system('mkdir '+path_direc+'/Equivalences')
	os.system('mkdir '+path_direc+'/OutPutFiles')
	os.system('mkdir '+path_direc+'/Frustration')
	#os.system('cp '+path_to_r+'/*.R* '+path_direc+'/Equivalences/')
	
def pdb_list(fasta_file):
	'''     This function is for the creation of folders and copying of necessary files
		Parameters:
			- fasta_file: MSA of input
			- pathPDB: path to the folder with the pdbs ()
	'''
	msa=open(fasta_file)
	sp=fasta_file.split('.')
	out=open(sp[0]+'.list','w')
	for lmsa in msa.readlines():
		if lmsa[0] == '>':
			out.write(lmsa[1:])
	out.close()
	msa.close()
	return(sp[0]+'.list')

def changes(JobID,MSA_File):
	'''     This function make changes in the MSA and make a List
		Parameters:
			- JobID: the job name
			- MSA_File: the MSA
			
	'''
	path_direc='FrustraEvo_'+JobID
	out=open(path_direc+'/MSA_Clean.fasta','w')
	out_list=open(path_direc+'/PDB_List.txt','w')
	pathAlign=MSA_File
	for seq_record in SeqIO.parse(pathAlign, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq
		out.write('>'+seqid+'\n'+str(seq)+'\n')
		
	out.close()
	out_list.close()

def FrustraPDB(list_pdbs,JobID,pathPDB='None'):
	'''     This function is for the frustration calculation
		Parameters:
			- list_pdbs: the list with the Pdbs 
			- JobID: the job name
			- pathPDB: path to the folder with the pdbs ()
			
	If you do not have the pdbs structures, the pipeline will download them from the Protein Data Bank database
	'''
	path_direc='FrustraEvo_'+JobID
	pdbs=open(list_pdbs,'r')
	frustdir=path_direc+'/Frustration/'
	for line in pdbs.readlines():
		line=line.rstrip('\n')
		if pathPDB == 'None':
			os.system('cd '+frustdir+'/;wget \'http://www.rcsb.org/pdb/files/'+line+'.pdb\' -O '+JobID+'/Frustration/'+line+'.pdb')
		else:
			#line=line.lower()
			os.system('cp '+pathPDB+'/'+line+'.pdb '+frustdir+'/'+line+'.pdb')
	pdbs.close()
	directory=os.getcwd()+'/'
	frustra=open(frustdir+'FrustraR.R','w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <- \''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'singleresidue\', ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
	
  
def checks(JobID):
	'''     This function checks the frustration calculations
	Parameters:
		- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	MSA=open(path_direc+'/MSA_Clean.fasta','r')
	out=open(path_direc+'/MSA_Chk.fasta','w')
	out_list=open(path_direc+'/PDB_ListChk.txt','w')	
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
	

def DeleteGaps(JobID):
	'''     This function find the reference protein using HHMSearch
		Parameters:
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	AClean=open(path_direc+'/MSA_Clean.fasta','w')
	Asearch=open(path_direc+'/MSASearch.fasta','w')
	pathAlign=path_direc+'/MSA_Chk.fasta'
	vector = list()
	c=0
	for seq_record in SeqIO.parse(pathAlign, 'fasta'):
		c = c + 1
		vector.append(seq_record.seq)
	porc=c*0.6 
	cgaps=0
	aux=''
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
	
	for seq_record in SeqIO.parse(pathAlign, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq
		pr=0
		AClean.write('>'+seqid+'\n')
		Asearch.write('>'+seqid+'\n')
		for i in range(0,len(vectorgaps)):
			if vectorgaps[i]==0:
				if pr==60:
					AClean.write('\n')
					pr=0
				AClean.write(seq[i])
				if seq[i] != '-':
					Asearch.write(seq[i])
				pr+=1
		AClean.write('\n')
		Asearch.write('\n')
	AClean.close()
	Asearch.close()
	os.system('hmmbuild --amino '+path_direc+'/family.hmm'+' '+path_direc+'/MSA_Clean.fasta')
	os.system('hmmsearch --domtblout '+path_direc+'/out_family'+' '+path_direc+'/family.hmm '+path_direc+'/MSASearch.fasta')

	hmmsal = open(path_direc+'/out_family','r')
	ini=0
	pdbidsal=''
	for line in hmmsal.readlines():
		if line[0] != '#':
			spline = line.split(' ')
			pdbidsal = spline[0]
			break
	hmmsal.close()
	return pdbidsal


def prepare_file(JobID,prot_ref):
	'''     This function prepare files for the calculations of the Equivalences
	Parameters:
		- JobID: the job name
		- prot_ref: the reference protein
	'''
	path_direc='FrustraEvo_'+JobID
	MSA=open(path_direc+'/MSA_Chk.fasta','r')
	out=open(path_direc+'/MSA_Chk_Ref.fasta','w')
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
	os.system('cat '+path_direc+'/MSA_Chk.fasta >> '+path_direc+'/MSA_Chk_Ref.fasta')


def MissingRes(JobID):
	'''     This function find Missing Residues
		Parameters:
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	hash = {'CYS' : 'C', 'ASP' : 'D', 'SER' : 'S','GLN' : 'Q','LYS' : 'K','ILE' : 'I','PRO' : 'P','THR' : 'T','PHE' : 'F','ASN' : 'N','GLY' : 'G','HIS' : 'H','LEU' : 'L','ARG' : 'R','TRP' : 'W','ALA' : 'A','VAL' : 'V','GLU' : 'E','TYR' : 'Y','MET' : 'M','MSE' : 'B'}
	MSA=open(path_direc+'/MSA_Chk_Ref.fasta','r')
	for line in MSA.readlines():
		line=line.rstrip('\n')
		if line[0] == '>':
			pdbch=line[1:]
			spline=pdbch.split('-')
			count=0
			srch=open(JobID+'/Frustration/'+spline[0]+'.pdb','r')
			for lsrch in srch.readlines():
				sp_sch=busca.split(' ')
				if len(spline) >= 2:
					if sp_sch[0] == 'REMARK' and sp_sch[1]=='465' and sp_sch[7] == spline[1]:
						chain.append(ch)
						pos=lsrch[21]+lsrch[22]+lsrch[23]+lsrch[24]+lsrch[25]
						com.append(pos.lstrip())
						iden.append(hashsp_sch[6])
						count=count+1
				else:
					if sp_sch[0] == 'REMARK' and sp_sch[1]=='465':
						com.append(sp_sch[10])
						chain.append(sp_sch[7])
						count=count+1
			srch.close()
			if count!=0:
				os.system('cp '+JobID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue '+JobID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue')
				frst_sr=open(JobID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue','r')
				frst_sr_out=open(JobID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','w')
				cnt=0
				spl=''
				cnt2=0
				lfrst_sr=frst_sr.readline()
				frst_sr_out.write('#Res ChainRes DensityRes AA NativeEnergy DecoyEnergy SDEnergy FrstIndex\n')
				for lfrst_sr in frst_sr.readlines():
					spl= lfrst_sr.split()
					if spl[0] == '#Res':
						frst_sr_out.write(lfrst_sr)
					else:
						if spl[0] > com[cnt]:
							while (int(spl[0]) - 1 >= int(com[cnt]) and count>cnt2):
								frst_sr_out.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n')
								cnt=cnt+1
								cnt2=cnt2+1
								if cnt == count:
									cnt=cnt-1
									break
							
							frst_sr_out.write(lfrst_sr)
						else:
							frst_sr_out.write(lfrst_sr)
				if spl[0]<com[cnt]:
					while cnt<=count-1:
						frst_sr_out.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n')
						cnt=cnt+1
				frst_sr.close()
				frst_sr_out.close()
	align.close()

def ChangeAlign():
	'''     This function make changes in the alignment in the case if you hace missing residues
	Parameters: 
		- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	os.system('cp '+path_direc+'/MSA_Chk_Ref.fasta '+path_direc+'/MSA_Chk_Ref_aux.fasta')
	pathAlign=path_direc+'/MSA_Chk_Ref_aux.fasta'

	out_msa=open(path_direc+'/MSA_Chk_Ref_aux_1.fasta','w')
	out_msa_ref=open(path_direc+'/MSA_Chk_Ref.fasta','w')

	for seq_record in SeqIO.parse(pathAlign, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq
		pdbid=seqid
		out_msa.write('>'+seqid+'\n')	
		out_msa_ref.write('>'+seqid+'\n')
		d=len(sseqid)
		frst_sl=open(JobID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r')
		c=0
		if(d>2):
			c=int(sseqid[d-2])
		frst_sl.readline()
		lfrst_sl = frst_sl.readline()
		sp_frst_sl = lfrst_sl.split(' ')
		i=int(sl[0])
		while int(i)<int(c):
			lfrst_sl = frst_sl.readline()
			i = i + 1
		w=0
		while True:
			frst_sl_pos =  lfrst_sl.split(' ') 
			if seq[w] == '-' or seq[w] == 'X' or seq[w] == 'B':
				out_msa.write('-')
				out_msa_ref.write('-')
			elif len(frst_sl_pos) >= 4 and frst_sl_pos[4] == 'Missing':
					out_msa.write('-')
					out_msa_ref.write('Z')
					l = frst_sl.readline()
			else:
				out_msa.write(seq[w])
				out_msa_ref.write(seq[w])
				lfrst_sl = frst_sl.readline()				
			w = w + 1		

			if w == len(seq):	
				break
		out_msa_ref.write('\n')
		out_msa.write('\n')
	out_msa.close()
	out_msa_ref.close()
	frst_sl.close()
	os.system('rm '+path_direc+'/MSA_Chk_Ref_aux_1.fasta')
	os.system('rm '+path_direc+'/MSA_Chk_Ref_aux.fasta')
	

def FinalAlign(JobID):
	'''     This function create the final files
	Parameters: 
		- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	pathAlign=path_direc+'/MSA_Chk_Ref.fasta'
	out_msa=open(path_direc+'/MSA_Final.fasta','w')
	out_pos=open(path_direc+'/Positions','w')

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
		frst_sr=''
		longline=0
		if c!=0:
			out_msa.write('>'+pdb+'\n')
			out_pos.write('>'+pdb+'\n')
			frst_sr=open(path_direc+'/Frustration/'+pdb+'.done/FrustrationData/'+pdb+'.pdb_singleresidue','r')
			lfrst_sr=frst_sr.readline()
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
			lfrst_sr = frst_sr.readline()
			if r==1 and len(pdbid)>2:
				splitres=lfrst_sr.split(' ')
				r= int(splitres[0])
				a = int(pdbid[2])-1
				while (r<a):
					lfrst_sr=frst_sr.readline()
					r = r + 1
			q=0
			for j in range (0,tam):
				if vector[j] == 0:
					if seq[j] != '-':
						lfrst_sr = frst_sr.readline()
				else:
					if seq[j] == 'Z':
						out_msa.write('-')
						longline+=1
						q = q + 1
					else:
						out_msa.write(str(seq[j]))
						longline+=1
						q = q + 1
					if seq[j] == 'Z':
						out_pos.write('Z ')
					elif seq[j] == '-':
						out_pos.write('G ')
					else:
						if len(splitres) < 2:
							splitres = lfrst_sr.split(' ')
						lfrst_sr = frst_sr.readline()
						out_pos.write(str(splitres[0])+' ')
						splitres = lfrst_sr.split(' ')
						r = r + 1
						
			out_msa.write('\n')
			out_pos.write('\n')
	frst_sr.close()


	MSA_Long=open(path_direc+'/Equivalences/long.txt','w')
	MSA_Long.write(str(longline)+'\n')
	MSA_Long.close()
	os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/long.txt')

	out_msa.close()
	out_pos.close()
	

def Equivalences(JobID):
	'''     This function create the final files
		Parameters: 
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	pos=open(path_direc+'/Positions','r')

	for line in pos.readlines():
		line=line[:-1]
		if line[0] == '>':
			pdbch=line[1:]
			out=open(path_direc+'/Equivalences/Equival_'+pdbch+'.txt','w')
			frst_sr=open(path_direc+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r')
			frst_sr_line=frst_sr.readline()
		else:
			ter=0
			sline=line.split(' ')
			splitres=frst_sr_line.split(' ')
			tam=len(sline)
			while(ter<tam - 1):
				ter+=1
				if sline[ter-1] == 'G' or sline[ter-1] == 'Z':
					out.write(str(ter)+'\tN/A\tN/A\tN/A\tN/A '+splitres[1]+'\n')
					
				else:	
					frst_sr_line=frst_sr.readline()
					frst_sr_line = frst_sr_line.rstrip('\n')
					splitres=frst_sr_line.split(' ')
					if frst_sr_line != '':
						if len(splitres) < 7 and splitres[4] !='Missing':
							break
						if int(splitres[0]) < int(sline[ter-1]):
							while True:
								frst_sr_line=frst_sr.readline()
								frst_sr_line = frst_sr_line.rstrip('\n')
								splitres=frst_sr_line.split(' ')
								if splitres[0] == sline[ter-1] or len(frst_sr_line)<1:
									break
						if len(splitres) == 6:
							ter-=1
						if splitres[0] == sline[ter-1] and len(splitres) > 7:
							if float(splitres[7]) > 0.55:
								out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MIN '+splitres[1]+'\n')
							elif float(splitres[7]) < -1:
								out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MAX '+splitres[1]+'\n')
							else:
								out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'NEU '+splitres[1]+'\n')
	out.close()
	frst_sr.close()

	pos.close()
	
	
def FastaMod(JobID):
	'''     This function create the file for the sequence logo: 
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	fasta=path_direc+'/MSA_Final.fasta'
	out=open(path_direc+'/Equivalences/Logo.fasta','w')

	for seq_record in SeqIO.parse(fasta, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq
		out.write(str(seq)+'\n')
	out.close()


def LogoCheck(JobID):
	'''     This function checks the logo 
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	os.system('cd '+path_direc+'/Equivalences/;cat *Equival* > AllEquivalB.txt')
	res=open(path_direc+'/Equivalences/AllEquivalB.txt','r')
	out=open(path_direc+'/Equivalences/AllEquival.txt','w')
	for lres in res:
		lres = lres[:-1]
		splres= lres.split('\t')
		tam=len(splres)
		if tam>4:
			out.write(str(lres)+'\n')

	res.close()
	out.close()


def plots_logo(JobID,prot_ref,path_to_r):
	'''     This function makes the plots for the frustration logo: 
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	checks(JobID)
	prepare_file(JobID,prot_ref)
	FinalAlign(JobID)
	Equivalences(JobID)
	FastaMod(JobID)
	LogoCheck(JobID)
	os.system('Rscript '+path_to_r+'/Logo.R --dir '+os.getcwd()+'/'+path_direc+'/Equivalences/')
	add_ref(JobID,prot_ref)
	os.system('Rscript '+path_to_r+'/Generator.R --dir '+os.getcwd()+'/'+path_direc+'/Equivalences/')
	os.system('mv '+path_direc+'/Equivalences/HistogramFrustration.png'+' '+path_direc+'/OutPutFiles/FrustrationLogo'+JobID+'.png')
	os.system('cp '+path_direc+'/MSA_Final.fasta'+' '+path_direc+'/OutPutFiles/MSA_'+JobID+'.fasta')
	os.system('mv '+path_direc+'/EvoFrustra-log'+' '+path_direc+'/OutPutFiles/Error_'+JobID+'.log')
	clean_files(JobID)
	VScript(JobID)

def add_ref(JobID,prot_ref):
	'''     This function generate the .pml files: 
			- JobID: the job name
			- prot_ref: the reference protein
	'''
	path_direc='FrustraEvo_'+JobID
	res=open(path_direc+'/Equivalences/Equival_'+prot_ref+'.txt','r')
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
	os.system('cp '+path_direc+'/Equivalences/CharactPosData'+' '+path_direc+'/OutPutFiles/IC_SingleRes_'+JobID)

def VScript(JobID):
	'''     This function generate the .pml files: 
			- JobID: the job name
	'''
	path_direc='FrustraEvo_'+JobID
	list_chk= open(path_direc+'/AuxFiles/PDB_ListChk.txt','r')
	for line in list_chk:
		line=line.rstrip('\n')
		pdbid=line
		sal = open(path_direc+'/Data/'+pdbid+'.done/FrustrationData/'+pdbid+'.pml','w')
		Equ = open(path_direc+'/Data/'+pdbid+'.done/Equival_'+pdbid+'.txt','r')
		ECon = open(path_direc+'/OutPutFiles/IC_SingleRes_'+JobID,'r')
		sal.write('load '+pdbid+'.pdb\nhide all\nshow cartoon, all\nbg_color white\ncolor black, all')
		EstCon=ECon.readline()
		for EstCon in ECon.readlines():
			Equi= Equ.readline()
			splitE= Equi.split()
			splitEC= EstCon.split()
			if splitEC[14] == 'MAX' and float(splitEC[13])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor red,resi '+splitE[1]+' and chain '+splitE[5])
					
			if splitEC[14] == 'MIN' and float(splitEC[13])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor green,resi '+splitE[1]+' and chain '+splitE[5])

			if splitEC[14] == 'NEU' and float(splitEC[13])>0.5:
				if splitEC[0] != splitE[0]:
					while int(splitEC[0]) == int(splitE[0]):
						EstCon=ECon.readline()
				if splitE[1] != 'N/A':
					sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor gray,resi '+splitE[1]+' and chain '+splitE[5])
		sal.close()
		Equ.close()
		ECon.close()
	list_chk.close()
	
def CMaps_Mutational(JobID,path_to_r,prot_ref):
	'''     This function generate the CMaps for mutational index: 
			- JobID: the job name
			- path_to_r: path to the R files
	'''
	directory=os.getcwd()+'/'
	path_direc='FrustraEvo_'+JobID
	frustdir='FrustraEvo_'+JobID+'/Frustration/'
	if os.path.isdir(path_direc+'/CMaps'):
		print('The Path: '+path_direc+'/CMaps exists')
	else:
		os.system('mkdir '+path_direc+'/CMaps')
		os.system('cp '+path_direc+'/PDB_ListChk.txt '+path_direc+'/CMaps/PDB_ListChk.txt')
		os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
		os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
	frustra=open(frustdir+'FrustraR.R','w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <- \''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'mutational\', ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
	MSA_Long=open(path_direc+'/Equivalences/long.txt','r')
	long=MSA_Long.readline()
	long=long[:-1]
	MSA_Long.close()
	dir_total=directory+path_direc
	os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Mutational.py '+str(long)+' '+dir_total)
	os.system('Rscript '+path_to_r+'/IC_conts_Mut.R --dir '+os.getcwd()+'/'+path_direc+'/CMaps/')
	add_ref_Cmaps(JobID,prot_ref,'Mut')
	os.system('cp '+path_direc+'/CMaps/IC_Mut.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JobID+'_Mut.png')

	
def CMaps_Configurational(JobID,path_to_r,prot_ref):
	'''     This function generate the CMaps for Configurational index: 
			- JobID: the job name
			- path_to_r: path to the R files
	'''
	directory=os.getcwd()+'/'
	path_direc='FrustraEvo_'+JobID
	frustdir='FrustraEvo_'+JobID+'/Frustration/'
	if os.path.isdir(path_direc+'/CMaps'):
		print('The Path: '+path_direc+'/CMaps exists')
	else:
		os.system('mkdir '+path_direc+'/CMaps')
		os.system('cp '+path_direc+'/PDB_ListChk.txt '+path_direc+'/CMaps/PDB_ListChk.txt')
		os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
		os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
	frustra=open(frustdir+'FrustraR.R','w')
	frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <- \''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'configurational\', ResultsDir = ResultsDir)\n')
	frustra.close()
	os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
	MSA_Long=open(path_direc+'/Equivalences/long.txt','r')
	long=MSA_Long.readline()
	long=long[:-1]
	MSA_Long.close()
	dir_total=directory+path_direc
	os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Conf.py '+str(long)+' '+dir_total)
	os.system('Rscript '+path_to_r+'/IC_conts_Conf.R --dir '+os.getcwd()+'/'+path_direc+'/CMaps/')
	os.system('cp '+path_direc+'/CMaps/IC_Conf.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JobID+'_Conf.png')
	add_ref_Cmaps(JobID,prot_ref,'Conf')
	
def add_ref_Cmaps(JobID,prot_ref,Mode):
	'''     This function generate the .pml files: 
			- JobID: the job name
			- prot_ref: the reference protein
	'''
	path_direc='FrustraEvo_'+JobID
	res=open(path_direc+'/Equivalences/Equival_'+prot_ref+'.txt','r')
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
	os.system('cp '+path_direc+'/CMaps/IC_'+Mode+'_ref'+' '+path_direc+'/OutPutFiles/IC_'+Mode+'_'+JobID)

def clean_files(JobID):
	path_direc='FrustraEvo_'+JobID
	os.system('mkdir '+path_direc+'/AuxFiles')
	os.system('rm '+path_direc+'/MSA_Final.fasta')
	os.system('mv '+path_direc+'/long.txt '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/MSA_Chk.fasta '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/MSA_Chk_Ref.fasta '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/MSA_Clean.fasta '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/PDB_List.txt '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/PDB_ListChk.txt '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/Positions '+path_direc+'/AuxFiles/')
	os.system('mv '+path_direc+'/Frustration/ '+path_direc+'/Data/')
	os.system('cd '+path_direc+'/Data;rm *.pdb*')
	lista=open(path_direc+'/AuxFiles/PDB_ListChk.txt')
	for line in lista.readlines():
		line=line.rstrip('\n')
		os.system('mv '+path_direc+'/Equivalences/Equival_'+line+'.txt '+path_direc+'/Data/'+line+'.done/')
	lista.close()
	os.system('rm -r '+path_direc+'/Equivalences/')
