# load required packages
suppressMessages(require(ggplot2))
suppressMessages(require(ggseqlogo))
suppressMessages(library(cowplot))
suppressMessages(library(seqinr))
suppressMessages(library(data.table))
suppressPackageStartupMessages(library("argparse"))  
parser <- ArgumentParser()
parser$add_argument("--dir", help="Path to CharactPosData file ")
args <- parser$parse_args()
FileTabla2=read.table(file=paste(args$dir, "CharactPosData", sep=""), header=T, stringsAsFactors=F)
msa = read.fasta(file=paste(args$dir, "MSA_Final.fasta", sep=""), forceDNAtolower = F, seqtype = "AA")

write.csv2(FileTabla2, file = paste(args$dir, "IC.csv", sep=""),row.names=T,col.names=T)  

h_minimos=as.numeric(FileTabla2[,11])
h_neutros=as.numeric(FileTabla2[,12])
h_maximos=as.numeric(FileTabla2[,13])

In=as.numeric(FileTabla2[1,3])
N=length(h_minimos) 

plot(h_minimos)
plot(h_maximos)
plot(h_neutros)

#generate data for bars
df_plot = data.table(N=rep(1:N,3),
                     h= c(h_minimos, h_maximos,h_neutros),
                     state= rep(c('Minimum', "Maximum", "Neutral"), each=N),
                     group=1)
start = 1
end = N

# generate logo
slog = ggseqlogo(unlist(getSequence(getFrag(msa, start,end), as.string = T)),
                 ncol = 100,
                 font= "akrobat_bold")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text())+
  scale_x_continuous(breaks = c(1,seq(5,2998, by = 5)), expand = c(0, 0))

# generate barplot
bplot = ggplot(subset(df_plot, group == 1 ), aes(x=N, y =h, fill = state))+
  geom_bar(stat= "identity",position="stack", color = "black")+
  scale_fill_manual(values= c("red4", "darkgreen", "grey"))+
  scale_x_continuous(breaks = c(1,seq(5,1998, by = 5)), expand = c(0, 0))+
  theme_logo()+
  xlab("")+
  ylab("")+
  theme(axis.text = element_text(size = 12),
        axis.ticks = element_line())

legend_logo = get_legend(slog)
legend_bplot = get_legend(bplot)
plot_grid(slog+theme(legend.position = "none"),
          bplot+theme(legend.position = "none"),
          ncol = 1, rel_heights = c(0.75,1), align= "v",axis = "lr")
# save figure
scale = end / 5
wd= N*20
png(paste(args$dir, "HistogramFrustration.png", sep=""),width = wd, height = 480, units = "px")
plot_grid(
  plot_grid(NULL,NULL,legend_logo, legend_bplot,NULL,NULL, ncol = 6) ,
  plot_grid(slog+theme(legend.position = "none"),
            bplot+theme(legend.position = "none"),
            ncol = 1, rel_heights = c(0.75,1), align= "v",axis = "lr"),
  rel_heights = c(0.1,1.5),
  ncol = 1)
dev.off()
