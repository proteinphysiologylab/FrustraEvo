res=open(res, sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/AllSalidaSResB.txt','r')
sal=open(res, sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/AllSalidaSRes.txt','w')

for lres in res:
	lres = lres[:-1]
	splres= lres.split(' ')
	tam=len(splres)
	if tam>5:
		sal.write(res+'\n')

res.close()
sal.close()
