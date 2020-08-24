import sys
import os

align=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/SeqAlign2.fasta','r')
sal=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/ListaPDBC.txt','w')

hash = {'CYS' : 'C', 'ASP' : 'D', 'SER' : 'S'
,'GLN' : 'Q','LYS' : 'K','ILE' : 'I','PRO' : 'P',
'THR' : 'T','PHE' : 'F','ASN' : 'N','GLY' : 'G','HIS' : 'H','LEU' : 'L','ARG' : 'R',
'TRP' : 'W','ALA' : 'A','VAL' : 'V','GLU' : 'E','TYR' : 'Y','MET' : 'M','MSE' : 'B'}

line=''
pdb=''
ch=''
pdbch=''
com=''
iden=[]
chain=[]
count=0
comi=''

for alin in align.readlines():
	alin=alin[:-1]
	cline=0
	lon=0
	com=[]
	iden=[]
	chain=[]
	sp=''
	ba=0
	ta=0
	if alin[0] == '>':
		line=alin[1:]
		spline=line.split('_')
		if len(spline) >= 2:
			sal.write(line)
			ch=spline[1]
			pdbch=spline[0]+'_'+spline[1]
		else:
			sal.write(spline[0]+'\n')
			pdbch=spline[0]
		splig= alin.split('_')
		tc=splig
		pdb=spline[0]
		count=0
		bus=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/'+pdb+'.pdb','r')
		for busca in bus.readlines():
			spbusca=busca.split(' ')
			if len(spline) >= 2:
				if spbusca[0] == 'REMARK' and spbusca[1]=='465' and spbusca[7] == ch:
					chain.append(ch)
					com.append(spbusca[10])
					iden.append(hash[spbusca[6]])
					count=count+1
			else:
				if spbusca[0] == 'REMARK' and spbusca[1]=='465':
					com.append(spbusca[10])
					chain.append(spbusca[7])
					count=count+1
		bus.close()
		if count==0:
			cp='cp '+sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdb+'.pdb_singleresidue '+sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdbch+'.pdb_msingleresidue'
			os.system(cp)
		else:
			sres=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdb+'.pdb_singleresidue','r')
			sressal=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdbch+'.pdb_msingleresidue','w')
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
sal.close()
