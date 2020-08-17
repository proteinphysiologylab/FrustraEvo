  FileTabla2=read.table(file="CharactPosDataN", header=T, stringsAsFactors=F)
  write.csv2(FileTabla2, file = "IC.csv",row.names=T,col.names=T)  
   
  h_minimos=as.numeric(FileTabla2[,8])
  h_neutros=as.numeric(FileTabla2[,9])
  h_maximos=as.numeric(FileTabla2[,10])
  
  N=length(h_minimos)
  
  plot(h_minimos)
  plot(h_maximos)
  plot(h_neutros)
  
  svg("HistogramFrustration.svg", height = 4, width =50)
  par(mar=c(2.5,2.5,2,2))
  par(mgp=c(2,0.5,0))
  barplot(names.arg=seq(from=1, to=N, by=1),  rbind(h_minimos, h_neutros, h_maximos), col=c("green", "gray", "red"), ylim=c(-0.30,2), axis.lty=1)
  box(lwd=2)
  legend(x="topleft", legend=c("Minimum", "Neutral", "Maximum"), pch=c(15, 15, 15), col=c("green", "gray", "red"))
  dev.off()

  