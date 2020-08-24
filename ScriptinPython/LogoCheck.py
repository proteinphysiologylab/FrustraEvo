import sys
import os

res=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/AllSalidaSResB.txt','r')
sal=open(sys.argv[1]+'/OutPutFiles'+sys.argv[2]+'/Equivalences/AllSalidaSRes.txt','w')

for lres in res:
	lres = lres[:-1]
	splres= lres.split('\t')
	tam=len(splres)
	if tam>4:
		sal.write(str(lres)+'\n')

res.close()
sal.close()
