  FileTabla2=read.table(file="CharactPosData", header=T, stringsAsFactors=F)
  write.csv2(FileTabla2, file = "IC.csv",row.names=T,col.names=T)  
   
  h_minimos=as.numeric(FileTabla2[,11])
  h_neutros=as.numeric(FileTabla2[,12])
  h_maximos=as.numeric(FileTabla2[,13])
  
  In=as.numeric(FileTabla2[1,3])
  N=length(h_minimos) + In - 1
  
  plot(h_minimos)
  plot(h_maximos)
  plot(h_neutros)
  
  svg("HistogramFrustration.svg", height = 4, width =50)
  par(mar=c(2.5,2.5,2,2))
  par(mgp=c(2,0.5,0))
  barplot(names.arg=seq(from=In, to=N, by=1),  rbind(h_minimos, h_neutros, h_maximos), col=c("darkgreen", "gray", "red4"), ylim=c(-0.30,2), axis.lty=1)
  box(lwd=2)
  legend(x="topleft", legend=c("Minimum", "Neutral", "Maximum"), pch=c(15, 15, 15), col=c("darkgreen", "gray", "red4"))
  dev.off()

  