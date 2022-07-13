tabla=read.table("AllSalidaSRes.txt", stringsAsFactors=F, header=F)
N=read.delim("long.txt", stringsAsFactors=F, header=F)


# limites frustracion single res
lim_max=-1
lim_min=0.55

frust_index=as.numeric(tabla[,4])#4 en el mio

hist(frust_index)

#m=matrix(nrow=N, ncol=N, data=0)

sink("CharactPosDataN")
cat("Res", "\t", "%Min", "\t" , "%Neu","\t","%Max","\t","CantMin","\t","CantNeu","\t","CantMax","\t","ICMin","\t","ICNeu","\t","ICMax","\t","ICTot","\t","FrustEstado","\n")




for (CharactPos in 1:N$V1)
{
  NOM=tabla[,2]
  minimos=length(which(tabla[,1]==CharactPos & tabla[,5]=='MIN'))  
  maximos=length(which(tabla[,1]==CharactPos & tabla[,5]=='MAX'))  
  neutros=length(which(tabla[,1]==CharactPos & tabla[,5]=='NEU'))
  
  minimos_t=minimos #total
  maximos_t=maximos
  neutros_t=neutros
  
  total=minimos+maximos+neutros
  
  #Laplace correction -> No es aplicable aca, porque se usa para predicciones. Si k (resultado)=0, entonces P= 1/n+2 donde n es el num de seqs usadas.
  #neutros[which(neutros==0)]=1/(total+2) 
  #maximos[which(maximos==0)]=1/(total+2) 
  #minimos[which(minimos==0)]=1/(total+2)
  
  
  #small sample correction: restar a Hbackground/esperada (s-1)/2*ln(2)*n donde n es el numero de secuencias sampleadas 
  correccion= (3-1)/(2*log(2)*total)
  
  minimos=minimos/total
  maximos=maximos/total
  neutros=neutros/total
  
  min_esperado=0.4 #frecuencias esperada
  max_esperado=0.1
  neu_esperado=0.5
  
  minimos_shannon=minimos*log2(minimos)
  maximos_shannon=maximos*log2(maximos)
  neutros_shannon=neutros*log2(neutros)
  
  #si no hay ningun residuo con ese estado de frustra, la entropia de ese estado es cero
  minimos_shannon[which(minimos==0)]=0
  maximos_shannon[which(maximos==0)]=0
  neutros_shannon[which(neutros==0)]=0
  
  shannon=-(minimos_shannon+neutros_shannon+maximos_shannon)
  h_shannon= -(0.4*log2(0.4)+0.1*log2(0.1)+0.5*log2(0.5)) -shannon
  h_shannon_corregida= h_shannon - correccion
  
  #   h_shannon= -(0.4*log2(0.4)+0.1*log2(0.1)+0.5*log2(0.5)) -shannon-((3-1)/(2*log(2)*total))
  
  #   h_shannon=log2(3)-shannon-((3-1)/(2*log(2)*total))
  
  
  h_minimos=minimos*h_shannon_corregida
  h_maximos=maximos*h_shannon_corregida
  h_neutros=neutros*h_shannon_corregida
  
  h_total=h_minimos+h_maximos+h_neutros
  
  if(maximos_t>neutros_t){
    if(maximos_t>minimos_t){estado='MAX'}
    else{estado='MIN'}
  }
  else{
    if(neutros_t>minimos_t){estado='NEU'}
    else{estado='MIN'}
  }
  
  wtf=CharactPos
  cat(wtf, "\t", minimos, "\t" , neutros,"\t",maximos,"\t",minimos_t,"\t",neutros_t,"\t",maximos_t,"\t",h_minimos,"\t",h_neutros,"\t",h_maximos,"\t",h_shannon_corregida,"\t",estado,"\n")
  
}


#sink()
