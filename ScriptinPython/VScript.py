
Lista= open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/ListaPDB.txt'

for lista in Lista:
	lista=lista[:-1]
	splitter=lista.split('_')
	sal = open(sys.argv[1]+'/OutPut'+sys.argv[2]+'/VisualizationScript/'+splitter[0]+'_'+splitter[1]+'.pml','w')
	Equ = open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/SalidaSRes'++splitter[0]+'_'+splitter[1]+'.txt','r')
	ECon = open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/CharactPosDataN','r')
	sal.write('load '+sys.argv[1]+'/OutPut'+sys.argv[2]+'/PDB/'+splitter[0]+'.pdb\nhide all\nshow cartoon, all\nbg_color white\ncolor black, all')
	EstCon=ECon.readline()
	for EstCon in ECon.readlines():
		Equi= Equ.readline()
		splitE= Equi.split(" ")
		splitEC= EstCon.split(" ")
		if splitEC[11] == 'MAX' and splitEC[10]>0.5:
			if splitEC[0] != splitE[0]:
				while splitEC[0] == splitE[0]:
					EstCon=ECon.readline()
			if splitE[1] != 'N/A':
				sal.write('\nshow sticks, resi '+splitE[1]+'\ncolor red,resi '+splitE[1])
				
		if splitEC[11] == 'MIN' and splitEC[10]>0.5:
			if splitEC[0] != @splitE[0]:
				while splitEC[0] == splitE[0]:
					EstCon=ECon.readline()
			if splitE[1] != 'N/A':
				sal.write('\nshow sticks, resi '+splitE[1]+'n\ncolor green,resi '+splitE[1]')

		if splitEC[11] == 'NEU' and splitEC[10]>0.5 :
			if splitEC[0] != @splitE[0]:
				while splitEC[0] == @splitE[0]:
					EstCon=ECon.readline()
			if splitE[1] != 'N/A':
				sal.write('\nshow sticks, resi '+splitE[1]+'n\ncolor gray,resi '+splitE[1]')
	sal.close()
	Equ.close()
	ECon.close()
Lista.close()
