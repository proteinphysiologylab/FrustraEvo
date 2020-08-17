
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
iden=''
chain=''
count=0;
comi=''

for alin in align.readlines():
	alin=align[:-1]
	cline=0
	lon=0
	com=''
	sp=''
	ba=0
	ta=0
	if alin[0] == '>':
		line=alin[:1]
		spline=line.split('_')
		sal.write(spline[0]+'_'+spline[1]+'\n'
		pdb=spline[0]
		ch=spline[1]
		pdbch=spline[0]+'_'+spline[1]
		splig= alin.split('_')
		comi=splig[2]
		tc=splig
		count=0
		bus=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+sp[0]+'.pdb','r');
		for busca in bus.readlines():
			spbusca=busca.split(' ')
			if spbusca[0] == 'REMARK' and spbusca[1]==465 and spbusca[3] == ch:
				com[count]=busca[4]
				chain[count]=busca[3]
				iden[count]=hash[busca[2]]
				count=count+1
		bus.close()
		if count==0:
			cp='cp '+sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdbch+'.pdb_singleresidue '+sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdb+'.pdb_singleresidue'
			os.system(cp)
		else:
			sres=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdbch+'.pdb_singleresidue','r');
			sressal=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Frustration/'+pdbch+'.pdb.done/FrustrationData/'+pdbch+'.pdb_singleresidue','w');
				cnt=0
				spl=''
				Sres=sres.readline()
				sressal.write(Sres)
				for Sresin sres.readlines():
					spl= Sres.split()
					if spl[0] > com[cnt]:
						while spl[0]-1 >= com[cnt] and count>cnt:
							sressal.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n';
							cnt=cnt+1
						sressal.write(Sres)
					else:
						sressal.write(Sres)
				if spl[0]<com[cnt]:
					while cnt<=count-1:
							sressal.write(com[cnt]+' '+chain[cnt]+' 0.000 '+iden[cnt]+' Missing Residue\n';
							cnt=cnt+1					
			sres.close()
			sressal.close()
align.close()
sal.close()
