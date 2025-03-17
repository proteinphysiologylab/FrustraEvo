#!/usr/bin/python3
import os
import re
import sys
from Bio import SeqIO
import numpy as np
#from scipy.cluster.hierarchy import linkage, dendrogram
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def pml_contactos(JodID,ref):
        path_direc='FrustraEvo_'+JodID
        mode=['configurational','mutational']
        mode2=['Conf','Mut']
        for i in range(0,len(mode)):
                list_chk = open(path_direc+'/AuxFiles/PDB_ListChk.txt','r')
                for line in list_chk.readlines():
                        pdbid=line.rstrip('\n')
                        if pdbid==ref:
                                table=open(path_direc+'/Data/Conservedresidues_'+mode[i],'w')
                        out = open(path_direc+'/Data/'+pdbid+'.done/VisualizationScrips/'+pdbid+'_contacts_'+mode[i]+'.pml','w')
                        out.write('load '+pdbid+'.pdb, '+pdbid.replace('-','_')+'\nhide line,'+pdbid.replace('-','_')+'\nunset dynamic_measures\nshow cartoon,'+pdbid.replace('-','_')+'\ncolor gray,'+pdbid.replace('-','_')+'\nrun draw_links.py\n')
                        IC = open(path_direc+'/OutPutFiles/IC_'+mode2[i]+'_'+JodID,'r')
                        for line in IC.readlines():
                        	line=line.rstrip('\n')
                        	sp=line.split()
                        	if sp[0] !='Res1':
	                        	if float(sp[21]) > 0.5 and float(sp[10]) >= 0.5:
                                                col=''
                                                if sp[22] == 'MAX':
                                                        col='red'
                                                elif sp[22] == 'MIN':
                                                        col='green'
                                                mut=open(path_direc+'/Data/'+pdbid+'.done/VisualizationScrips/'+pdbid+'.pdb_'+mode[i]+'.pml','r')
                                                for lmut in mut.readlines():
                                                        lmut=lmut.rstrip('\n')
                                                        spmut=lmut.split()
                                                        if sp[4] in lmut and sp[6] in lmut:
                                                                if pdbid==ref and col !='':
                                                                      table.write(sp[4]+' '+sp[6]+' '+col+'\n')
                                                                if 'green' in lmut and col == 'red':
                                                                      lmut=lmut.replace('green',col)
                                                                      lmut=lmut.replace('min_frst','max_frst')
                                                                elif 'red' in lmut and col == 'green':
                                                                      lmut=lmut.replace('red',col)
                                                                      lmut=lmut.replace('max_frst','min_frst')
                                                                out.write(lmut+'\n')
                                                                break
                                                mut.close()
                        if pdbid==ref:
                                table.close() 
                        out.write('zoom all\nhide labels\ncolor red, max_frst_wm_'+pdbid.replace('-','_')+'\ncolor green, min_frst_wm_'+pdbid.replace('-','_'))

                        IC.close()
                        out.close()
                        list_chk.close()

def copyfiles(JodID,path_to_r,path_to_Pdbs):

        '''     This function is for the creation of folders and copying of necessary files
                Parameters:
                        - JodID: the job name
                        - path_to_r: path to the R files
        '''
        path_direc='FrustraEvo_'+JodID
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

def changes(JodID,MSA_File):
        '''     This function make changes in the MSA and make a List
                Parameters:
                        - JodID: the job name
                        - MSA_File: the MSA

        '''
        path_direc='FrustraEvo_'+JodID
        out=open(path_direc+'/MSA_Clean_aux.fasta','w')
        out_list=open(path_direc+'/PDB_List.txt','w')
        pathAlign=MSA_File
        for seq_record in SeqIO.parse(pathAlign, 'fasta'):
                seqid=seq_record.id
                seq=seq_record.seq
                out.write('>'+seqid+'\n'+str(seq.upper()).replace('X','-')+'\n')

        out.close()
        out_list.close()

def FrustraPDB(list_pdbs, JodID, seqdist, pathPDB='None'):
        '''     This function is for the frustration calculation
                Parameters:
                        - list_pdbs: the list with the Pdbs 
                        - JodID: the job name
                        - pathPDB: path to the folder with the pdbs ()

        If you do not have the pdbs structures, the pipeline will download them from the Protein Data Bank database
        '''
        path_direc='FrustraEvo_'+JodID
        frustdir=path_direc+'/Frustration/'
        directory=os.getcwd()+'/'
        frustra=open(frustdir+'FrustraR.R','w')
        #frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <-\''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'singleresidue\', ResultsDir = ResultsDir)\n')
        #added by mfernan3 to include the seqdist flag
        frustra.write(f"library(frustratometeR)\n"
              f"PdbsDir <- '{directory}{frustdir}'\n"
              f"ResultsDir <- '{directory}{frustdir}'\n"
              f"dir_frustration(\n"
              f"    PdbsDir = PdbsDir,\n"
              f"    SeqDist = {seqdist},\n"
              f"    Mode = 'singleresidue',\n"
              f"    ResultsDir = ResultsDir\n"
              f")\n")
        frustra.close()
        #os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
        #added by mfernan3 to supress PDBConstructionWarnings
        command = f'cd {frustdir} && Rscript FrustraR.R > FrustraR.log 2>&1'
        os.system(command)


def obtain_seq(pdbid,JodID):
	dic =  {'ALA' : 'A','ARG' : 'R','ASN' : 'N','ASP' : 'D','CYS' : 'C','GLN' : 'Q','GLU' : 'E','GLY' : 'G','HIS' : 'H','ILE' : 'I','LEU' : 'L','LYS' : 'K','MET' : 'M','PHE' : 'F','PRO' : 'P','SER' : 'S','THR' : 'T','TRP' : 'W','TYR' : 'Y', 'VAL' : 'V'}
	seq=''
	path_direc='FrustraEvo_'+JodID
	pdb=open(path_direc+'/Frustration/'+pdbid+'.pdb')
	nn=-100
	for lpdb in pdb.readlines():
		if 'ATOM' == lpdb[0:4] and len(lpdb) > 60:
			aa=lpdb[17]+lpdb[18]+lpdb[19]
			if nn != int(lpdb[22:26]) and lpdb[16]==' ' and lpdb[26]==' ' and aa in dic:
				aa=lpdb[17]+lpdb[18]+lpdb[19]
				seq+=dic[aa]
			nn = int(lpdb[22:26])
	pdb.close()
	return seq

def checks_seq(list_pdbs,JodID,pathPDB='None'):
	path_direc='FrustraEvo_'+JodID
	seqs=0
	frustdir=path_direc+'/Frustration/'
	MSA=open(path_direc+'/MSA_Clean_aux.fasta','r')
	out=open(path_direc+'/MSA_Clean.fasta','w')
	#out_log=open(path_direc+'/ErrorSeq.log','w')
        #len_seq=0
	for seq_record in SeqIO.parse(MSA, 'fasta'):
		seqid=seq_record.id
		seq=seq_record.seq.replace('-','')
		seq_to_print=seq_record.seq
		pathpdb=path_direc+'/Frustration/'+seqid+'.pdb'
		pathdb=pathPDB+'/'+seqid+'.pdb'
		if not os.path.exists(pathdb):
			os.system('cd '+frustdir+'/;wget \'http://www.rcsb.org/pdb/files/'+seqid+'.pdb\' -O '+JodID+'/Frustration/'+seqid+'.pdb')
		else:
			os.system('cp '+pathPDB+'/'+seqid+'.pdb '+frustdir+'/'+seqid+'.pdb')
		if os.path.exists(pathpdb):
			seqpdb=obtain_seq(seqid,JodID)
			if seq == seqpdb:
				print('Sequence '+seqid+' checked')
				out.write('>'+str(seqid)+'\n'+str(seq_to_print)+'\n')
			else:
				out_log=open(path_direc+'/ErrorSeq.log','a')
				if seqs == 0:
					out_log.write('There are errors in the input data, fix them and submit the job again. Take into account that IDs are matched in a case-sensitive way. Take into account that the sequences need to be aligned.\n')
				print(seqid+' sequence does not match between the MSA and the pdb file')
				#out_log=open(path_direc+'/ErrorSeq.log','a')
				out_log.write(seqid+' sequence does not match between the MSA and the pdb file\n')
				os.system('rm '+frustdir+'/'+seqid+'.pdb')
				seqs+=1
				out_log.close()
	MSA.close()
	out.close()
#	archivo_fasta = path_direc+'/MSA_Clean_aux.fasta'
#	msa_c, n_seq = es_MSA_valido(archivo_fasta)
#	if not msa_c:
#		if n_seq == 1:
#			out_log=open(path_direc+'/ErrorSeq.log','a')
#			out_log.write('The Fasta file provided is not a valid fasta format.\n')
#		else:
#			out_log=open(path_direc+'/ErrorSeq.log','a')
#			out_log.write('The Fasta file provided has less than 3 sequences, the minimum number of sequences is 3.\n')
	#out_log.write(str(seqs))
	#out_log.close()ç
	#if seqs == 0:
	#	os.sytem('rm '+path_direc+'/ErrorSeq.log')
	#return seqs

def es_MSA_valido(archivo_fasta):
    # Cargar el archivo FASTA
	secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))

    # Verificar si todas las secuencias tienen la misma longitud
	c = 1
	longitud_referencia = len(secuencias[0])
	for seq in secuencias[1:]:
		c += 1
		if len(seq) != longitud_referencia:
			return False, 1
	if c < 3:
		return False, 0
	return True, 1

def checks(JodID):
        '''     This function checks the frustration calculations
        Parameters:
                - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
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


def DeleteGaps(JodID):
        '''     This function find the reference protein using HHMSearch
                Parameters:
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
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


def prepare_file(JodID,prot_ref):
        '''     This function prepare files for the calculations of the Equivalences
        Parameters:
                - JodID: the job name
                - prot_ref: the reference protein
        '''
        path_direc='FrustraEvo_'+JodID
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


def MissingRes(JodID):
        '''     This function find Missing Residues
                Parameters:
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        hash = {'CYS' : 'C', 'ASP' : 'D', 'SER' : 'S','GLN' : 'Q','LYS' : 'K','ILE' : 'I','PRO' : 'P','THR' : 'T','PHE' : 'F','ASN' : 'N','GLY' : 'G','HIS' : 'H','LEU' : 'L','ARG' : 'R','TRP' : 'W','ALA' : 'A','VAL' : 'V','GLU' : 'E','TYR' : 'Y','MET' : 'M','MSE' : 'B'}

        MSA=open(path_direc+'/MSA_Chk_Ref.fasta','r')
        for line in MSA.readlines():
                line=line.rstrip('\n')
                if line[0] == '>':
                        pdbch=line[1:]
                        spline=pdbch.split('-')
                        count=0
                        srch=open(JodID+'/Frustration/'+spline[0]+'.pdb','r')
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
                                os.system('cp '+JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue '+JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue')
                                frst_sr=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_msingleresidue','r')
                                frst_sr_out=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','w')
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
                - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
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
                frst_sl=open(JodID+'/Frustration/'+pdbch+'.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r')
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
	

def FinalAlign(JodID):
        '''     This function create the final files
        Parameters: 
                - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        pathAlign=path_direc+'/MSA_Chk_Ref.fasta'
        out_msa=open(path_direc+'/MSA_Final.fasta','w')
        out_pos=open(path_direc+'/Positions','w')

        vector = list()
        c=0     
        m=0
        control=0 
        cgaps=0
        tam=0
        tam=0
        r=0

        for seq_record in SeqIO.parse(pathAlign, 'fasta'):
                seqid=seq_record.id
                seq=seq_record.seq
                pdb=seqid
                longline=0
                if c!=0:
                        out_msa.write('>'+pdb+'\n')
                        out_pos.write('>'+pdb+'\n')
                        frst_sr=open(path_direc+'/Frustration/'+pdb+'.done/FrustrationData/'+pdb+'.pdb_singleresidue','r')
                        lfrst_sr=frst_sr.readlines()
                if c==0:
                        c=1
                        for i in range(0,len(seq)):
                                if seq[i] == '-' or seq[i] == 'Z':
                                        vector.append(0)
                                else:
                                        vector.append(1)
                else:
                        splitres=''
             #           if r==1 and len(pdbid)>2:
              #                  splitres=lfrst_sr.split(' ')
               #                 r= int(splitres[0])
                #                a = int(pdbid[2])-1
                 #               while (r<a):
                  #                      lfrst_sr=frst_sr.readline()
                   #                     r = r + 1
                        q=0
                        for j in range (0,len(seq)):
                                if vector[j] == 0:
                                        if seq[j] != '-':
                                                q += 1
                                else:
                                        if seq[j] == 'Z' or seq[j] == 'X':
                                                out_msa.write('-')
                                                longline+=1
                                                q += 1
                                        else:
                                                out_msa.write(str(seq[j]))
                                                longline+=1
                                                if seq[j] != '-':
                                                     q += 1
                                        if seq[j] == 'Z':
                                                out_pos.write('Z ')
                                        elif seq[j] == '-':
                                                out_pos.write('G ')
                                        else:
                                                splitres = lfrst_sr[q].split(' ')
                                                out_pos.write(str(splitres[0])+' ')
                                                r += 1

                        out_msa.write('\n')
                        out_pos.write('\n')
                        frst_sr.close()
        MSA_Long=open(path_direc+'/Equivalences/long.txt','w')
        MSA_Long.write(str(longline)+'\n')
        MSA_Long.close()
        os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/long.txt')

        out_msa.close()
        out_pos.close()


def Equivalences(JodID):
        '''     This function create the final files
                Parameters: 
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        pos=open(path_direc+'/Positions','r')

        for line in pos.readlines():
                line=line.rstrip('\n')
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
                                        out.write(str(ter)+'\tN/A\tN/A\tN/A\tN/A\t'+splitres[1]+'\n')

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
                                                                out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MIN\t'+splitres[1]+'\n')
                                                        elif float(splitres[7]) < -1:
                                                                out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'MAX\t'+splitres[1]+'\n')
                                                        else:
                                                                out.write(str(ter)+'\t'+str(splitres[0])+'\t'+str(splitres[3])+'\t'+str(splitres[7])+'\t'+'NEU\t'+splitres[1]+'\n')
        out.close()
        frst_sr.close()

        pos.close()


def FastaMod(JodID):
        '''     This function create the file for the sequence logo: 
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        fasta=path_direc+'/MSA_Final.fasta'
        out=open(path_direc+'/Equivalences/Logo.fasta','w')

        for seq_record in SeqIO.parse(fasta, 'fasta'):
                seqid=seq_record.id
                seq=seq_record.seq
                out.write(str(seq)+'\n')
        out.close()


def LogoCheck(JodID):
        '''     This function checks the logo 
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
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


def plots_logo(JodID,prot_ref,path_to_r):
        '''     This function makes the plots for the frustration logo: 
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        os.system('cp '+path_direc+'/MSA_Final.fasta'+' '+path_direc+'/Equivalences/MSA_Final.fasta')
        os.system('Rscript '+path_to_r+'/Logo.R --dir '+os.getcwd()+'/'+path_direc+'/Equivalences/')
        add_ref(JodID,prot_ref)
        os.system('Rscript '+path_to_r+'/Generator.R --dir '+os.getcwd()+'/'+path_direc+'/Equivalences/')
        os.system('mv '+path_direc+'/Equivalences/HistogramFrustration.png'+' '+path_direc+'/OutPutFiles/FrustrationLogo'+JodID+'.png')
        os.system('cp '+path_direc+'/MSA_Final.fasta'+' '+path_direc+'/OutPutFiles/MSA_'+JodID+'.fasta')
        os.system('mv '+path_direc+'/EvoFrustra-log'+' '+path_direc+'/OutPutFiles/Error_'+JodID+'.log')
        #os.system('tar -zcvf '+JodID+'.tar.gz '+path_direc)

def add_ref(JodID,prot_ref):
        '''     This function generate the .pml files: 
                        - JodID: the job name
                        - prot_ref: the reference protein
        '''
        path_direc='FrustraEvo_'+JodID
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
        out.write('Res\tAA_Ref\tNum_Ref\tProt_Ref\t%Min\t%Neu\t%Max\tCantMin\tCantNeu\tCantMax\tICMin\tICNeu\tICMax\tICTot\tFrustIC\n')
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

def VScript(JodID,tiempo):
        '''     This function generate the .pml files: 
                        - JodID: the job name
        '''
        path_direc='FrustraEvo_'+JodID
        out=open(path_direc+'/Time.log','w')
        out.write(str(tiempo)+'\n')
        out.close()
        list_chk= open(path_direc+'/AuxFiles/PDB_ListChk.txt','r')
        for line in list_chk:
                line=line.rstrip('\n')
                pdbid=line
                sal = open(path_direc+'/Data/'+pdbid+'.done/FrustrationData/'+pdbid+'.pml','w')
                Equ = open(path_direc+'/Data/'+pdbid+'.done/Equival_'+pdbid+'.txt','r')
                ECon = open(path_direc+'/OutPutFiles/IC_SingleRes_'+JodID,'r')
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

def CMaps_Mutational(JodID, seqdist, path_to_r,prot_ref):
        print('SeqdistMUT '+seqdist)
        '''     This function generate the CMaps for mutational index: 
                        - JodID: the job name
                        - path_to_r: path to the R files
        '''
        directory=os.getcwd()+'/'
        path_direc='FrustraEvo_'+JodID
        frustdir='FrustraEvo_'+JodID+'/Frustration/'
        if os.path.isdir(path_direc+'/CMaps'):
                print('The Path: '+path_direc+'/CMaps exists')
        else:
                os.system('mkdir '+path_direc+'/CMaps')
                os.system('cp '+path_direc+'/PDB_ListChk.txt '+path_direc+'/CMaps/PDB_ListChk.txt')
                os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
                os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
        frustra=open(frustdir+'FrustraR.R','w')
        #frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <- \''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'mutational\', ResultsDir = ResultsDir, Graphics = FALSE)\n')
        #added by mfernan3 to include the seqdist flag
        frustra.write(f"library(frustratometeR)\n"
              f"PdbsDir <- '{directory}{frustdir}'\n"
              f"ResultsDir <- '{directory}{frustdir}'\n"
              f"dir_frustration(\n"
              f"    PdbsDir = PdbsDir,\n"
              f"    SeqDist = {seqdist},\n"
              f"    Mode = 'mutational',\n"
              f"    ResultsDir = ResultsDir,\n"
              f"    Graphics = FALSE\n"
              f")\n")
        frustra.close()
        #os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
        #added by mfernan3 to supress PDBConstructionWarnings
        command = f'cd {frustdir} && Rscript FrustraR.R > FrustraR.log 2>&1'
        os.system(command)
        MSA_Long=open(path_direc+'/Equivalences/long.txt','r')
        long=MSA_Long.readline()
        long=long[:-1]
        MSA_Long.close()
        dir_total=directory+path_direc
        os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Mutational.py '+str(long)+' '+dir_total+' '+prot_ref)
        os.system('Rscript '+path_to_r+'/IC_conts_Mut.r --dir '+os.getcwd()+'/'+path_direc+'/CMaps/')
        add_ref_Cmaps(JodID,prot_ref,'Mut')
        os.system('cp '+path_direc+'/CMaps/IC_Mut.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JodID+'_Mut.png')


def CMaps_Configurational(JodID, seqdist, path_to_r,prot_ref):
        print('SeqdistCONF '+seqdist)
        '''     This function generate the CMaps for Configurational index: 
                        - JodID: the job name11
                        - path_to_r: path to the R files
        '''
        directory=os.getcwd()+'/'
        path_direc='FrustraEvo_'+JodID
        frustdir='FrustraEvo_'+JodID+'/Frustration/'
        if os.path.isdir(path_direc+'/CMaps'):
                print('The Path: '+path_direc+'/CMaps exists')
        else:
                os.system('mkdir '+path_direc+'/CMaps')
                os.system('cp '+path_direc+'/PDB_ListChk.txt '+path_direc+'/CMaps/PDB_ListChk.txt')
                os.system('cp '+path_direc+'/Equivalences/long.txt '+path_direc+'/CMaps/long.txt')
                os.system('cp '+path_to_r+'/*.py* '+path_direc+'/CMaps/')
        frustra=open(frustdir+'FrustraR.R','w')
        #frustra.write('library(frustratometeR)\nPdbsDir <- \''+directory+frustdir+'\'\nResultsDir <- \''+directory+frustdir+'\'\ndir_frustration(PdbsDir = PdbsDir, Mode = \'configurational\', ResultsDir = ResultsDir,Graphics = FALSE)\n')
        #added by mfernan3 to include the seqdist flag
        frustra.write(f"library(frustratometeR)\n"
              f"PdbsDir <- '{directory}{frustdir}'\n"
              f"ResultsDir <- '{directory}{frustdir}'\n"
              f"dir_frustration(\n"
              f"    PdbsDir = PdbsDir,\n"
              f"    SeqDist = {seqdist},\n"
              f"    Mode = 'configurational',\n"
              f"    ResultsDir = ResultsDir,\n"
              f"    Graphics = FALSE\n"
              f")\n")
              
        frustra.close()
        #os.system('cd '+frustdir+';Rscript FrustraR.R > FrustraR.log')
        #added by mfernan3 to supress PDBConstructionWarnings
        command = f'cd {frustdir} && Rscript FrustraR.R > FrustraR.log 2>&1'
        os.system(command)
        MSA_Long=open(path_direc+'/Equivalences/long.txt','r')
        long=MSA_Long.readline()
        long=long[:-1]
        MSA_Long.close()
        dir_total=directory+path_direc
        os.system('cd '+path_direc+'/CMaps;python3 IC_Conts_Conf.py '+str(long)+' '+dir_total+' '+prot_ref)
        os.system('Rscript '+path_to_r+'/IC_conts_Conf.r --dir '+os.getcwd()+'/'+path_direc+'/CMaps/')
        os.system('cp '+path_direc+'/CMaps/IC_Conf.png'+' '+path_direc+'/OutPutFiles/CMaps'+'_'+JodID+'_Conf.png')
        add_ref_Cmaps(JodID,prot_ref,'Conf')

def add_ref_Cmaps(JodID,prot_ref,Mode):
        '''     This function generate the .pml files: 
                        - JodID: the job name
                        - prot_ref: the reference protein
        '''
        path_direc='FrustraEvo_'+JodID
        res=open(path_direc+'/Equivalences/Equival_'+prot_ref+'.txt','r')
        out=open(path_direc+'/CMaps/IC_'+Mode+'_ref','w')
        ic=open(path_direc+'/CMaps/IC_'+Mode,'r')
        vectorAA=[]
        vectorchain=[]
        num=[]
        for line in res.readlines():
                line=line.rstrip('\n')
                sp=line.split()
                vectorAA.append(sp[2])
                vectorchain.append(sp[5])
                num.append(sp[1])
        out.write('Res1\tRes2\tAA1\tAA2\tNumRes1_Ref\tChain1_Ref\tNumRes2_Ref\tChain2_Ref\tProt_Ref\tNoContacts\tFreqConts\tpNEU\tpMIN\tpMAX\tHNEU\tHMIN\tHMAX\tHtotal\tICNEU\tICMIN\tICMAX\tICtotal\tFstConserved\n')
        a=0
        for line in ic.readlines():
                line=line.rstrip('\n')
                if a==0:
                        a+=1
                else:
                        sp=line.split()
                        out.write(sp[0]+'\t'+sp[1]+'\t'+vectorAA[int(sp[0])-1]+'\t'+vectorAA[int(sp[1])-1]+'\t'+num[int(sp[0])-1]+'\t'+vectorchain[int(sp[1])-1]+'\t'+str(num[int(sp[1])-1])+'\t'+str(vectorchain[int(sp[1])-1])+'\t'+prot_ref+'\t'+sp[2]+'\t'+sp[3]+'\t'+sp[4]+'\t'+sp[5]+'\t'+sp[6]+'\t'+sp[7]+'\t'+sp[8]+'\t'+sp[9]+'\t'+sp[10]+'\t'+sp[11]+'\t'+sp[12]+'\t'+sp[13]+'\t'+sp[14]+'\t'+sp[15]+'\n')

        ic.close()
        out.close()
        res.close()
        os.system('cp '+path_direc+'/CMaps/IC_'+Mode+'_ref'+' '+path_direc+'/OutPutFiles/IC_'+Mode+'_'+JodID)
        
def df_MSFA(JodID):
        path_direc='FrustraEvo_'+JodID
        pos=open(path_direc+'/AuxFiles/PDB_ListChk.txt','r')
        out=open(path_direc+'/AuxFiles/DF_Colores','w')
        out.write('PDBID,pos,AA,FstState,FstI\n')
        for line in pos.readlines():
                line=line.rstrip('\n')
                equival=open(path_direc+'/Data/'+line+'.done/Equival_'+line+'.txt','r')
                for lequiv in equival.readlines():
                        sp=lequiv.split('\t')
                        if sp[2] == 'N/A':
                               out.write(line+','+sp[0]+',-,-,N/A\n')
                        else:
                               out.write(line+','+sp[0]+','+sp[2]+','+sp[4]+','+sp[3]+'\n')
                equival.close()
        out.close()
        pos.close()


def clean_files(JodID,path_r,prot_ref, cmaps):
        path_direc='FrustraEvo_'+JodID
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
        #os.system('cd '+path_direc+'/Data;rm *.pdb*')
        if cmaps.lower() == 'yes': #add clause for cmaps
                os.system('rm -r '+path_direc+'/CMaps')
        #os.system('rm '+path_direc+'/ *.py*')
        lista=open(path_direc+'/AuxFiles/PDB_ListChk.txt')
        for line in lista.readlines():
                line=line.rstrip('\n')
                os.system('mv '+path_direc+'/Equivalences/Equival_'+line+'.txt '+path_direc+'/Data/'+line+'.done/')
        lista.close()
        os.system('rm -r '+path_direc+'/Equivalences/')
        #ACA LLAMAR A LA FUNCION
        df_MSFA(JodID)
        os.system('Rscript '+path_r+'/MSFA.R --dir '+os.getcwd()+'/'+path_direc+'/ --jobid '+JodID)
        if cmaps.lower() == 'yes': #add clause for cmaps
                pml_contactos(JodID,prot_ref)


