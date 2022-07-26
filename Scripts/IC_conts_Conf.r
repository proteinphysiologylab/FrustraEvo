suppressPackageStartupMessages(library("argparse"))  
parser <- ArgumentParser()
parser$add_argument("--dir", help="long.txt")
args <- parser$parse_args()
N=read.delim(paste(args$dir, "long.txt", sep=""), stringsAsFactors=F, header=F)
rangoresiduos<- c(1,N$V1)

IC_Conts_Conf <- read.delim(paste(args$dir, "IC_Conf", sep=""))

#------------------Configurational
maxfrust_Conf<- subset(IC_Conts_Conf, IC_Conts_Conf$EstadoConservado=='MAX')
minfrust_Conf<-subset(IC_Conts_Conf, IC_Conts_Conf$EstadoConservado=='MIN')
neufrust_Conf<-subset(IC_Conts_Conf, IC_Conts_Conf$EstadoConservado=='NEU')

max<-hist(maxfrust_Conf$ICtotal*maxfrust_Conf$FreqConts, breaks=20, probability=1)
min<-hist(minfrust_Conf$ICtotal*minfrust_Conf$FreqConts, breaks=20, probability=1)
neu<-hist(neufrust_Conf$ICtotal*neufrust_Conf$FreqConts, breaks=20, probability=1)


png(filename=paste(args$dir, "histICConf.png", sep=""), height=400, width=600, bg="white")
plot(neu, col='grey', xlab='IC of Confational Frustration')
lines(min, col='green')
lines(max, col='red')
dev.off()

densmax<-density(maxfrust_Conf$ICtotal*maxfrust_Conf$FreqConts, bw=0.07)
densmin<-density(minfrust_Conf$ICtotal*minfrust_Conf$FreqConts,bw=0.07)
densneu<-density(neufrust_Conf$ICtotal*neufrust_Conf$FreqConts,bw=0.07)

png(filename=paste(args$dir, "densICConf.png", sep=""), height=400, width=600, bg="white")
plot(densmax, col='red', xlab='FIC')
lines(densneu, col='grey')
lines(densmin, col='green')
dev.off()

#--------mapa de freq pintado en escala de grises- Conf
rbPal <- colorRampPalette(c('grey','black'))

## 10 = numero de grupos en los que se corta la frecuencia.
Col <- rbPal(10)[as.numeric(cut(IC_Conts_Conf$FreqConts,breaks = 10))]
cuts<-levels(cut(IC_Conts_Conf$FreqConts,breaks = 10))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","[",cuts)

png(filename=paste(args$dir, "IC_Conf.png", sep=""), height=1200, width=1200, bg="white")

# Add extra space to right of plot area; change clipping to figure
par(mar=c(6, 6, 4.1, 8), xpd=TRUE)

plot(IC_Conts_Conf$Res.1,IC_Conts_Conf$Res,pch=15, cex=1, col=Col, cex.axis=2, cex.lab=2, xlim=rangoresiduos, ylim=rangoresiduos, ylab="Res", xlab="Res")

points(neufrust_Conf$Res, neufrust_Conf$Res.1, col='grey', pch=15, cex=1)
points(minfrust_Conf$Res, minfrust_Conf$Res.1, col='green', pch=15, cex=1)
points(maxfrust_Conf$Res, maxfrust_Conf$Res.1, col='red', pch=15, cex=1)
lines(x=rangoresiduos, y=rangoresiduos, lty=2)

# Add legend to top right, outside plot region
legend("topright",inset=c(-0.22,0),cuts,col=rbPal(10),pch=16, cex=1)

dev.off()
