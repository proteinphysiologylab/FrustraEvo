suppressMessages(require(ggplot2))
suppressMessages(require(ggseqlogo))
suppressMessages(library(cowplot))
suppressMessages(library(seqinr))
suppressMessages(library(data.table))

msa = read.fasta(file = "./SalidaAlignSE.fasta", forceDNAtolower = F, seqtype = "AA")

FileTabla2=read.table(file="./CharactPosDataN", header=T, stringsAsFactors=F)

h_minimos=as.numeric(FileTabla2[,8])
h_neutros=as.numeric(FileTabla2[,9])
h_maximos=as.numeric(FileTabla2[,10])

N=length(h_minimos)


if (N>100){
  df_plot = data.table(N=rep(1:N,3),
                       h= c(h_minimos, h_maximos,h_neutros),
                       state= rep(c('Minimum', "Maximum", "Neutral"), each=N),
                       group=c(rep(1, 99), rep(2:100,each = 100,
                                               length.out = ifelse(N> 99, N-99, N))))
  start = c(1, seq(100,N, 100))
  end = c(seq(99,N, 100), N)
} else{
  df_plot = data.table(N=rep(1:N,3),
                       h= c(h_minimos, h_maximos,h_neutros),
                       state= rep(c('Minimum', "Maximum", "Neutral"), each=N),
                       group=1)
  start = 1
  end = N
}

list_plots = list()


for (i in 1:max(df_plot$group)){
  print (i)
  
  slog = ggseqlogo(unlist(getSequence(getFrag(msa, start[i],end[i]), as.string = T)),
                   ncol = 100,
                   font= "akrobat_bold")+
    theme(axis.text.x = element_blank())+
    ylab("SeqIC")
  bplot = ggplot(subset(df_plot, group == i ), aes(x=N, y =h, fill = state))+
    geom_bar(stat= "identity",position="stack", color = "black")+
    scale_fill_manual(values= c("red4", "darkgreen", "grey"))+
    scale_x_continuous(breaks = c(1,seq(10,1998, by = 10)))+
    theme_logo()+
    xlab("")+
    ylab("FrustIC")+
    theme(axis.text = element_text(size = 12))
  
  legend_logo = get_legend(slog)
  legend_bplot = get_legend(bplot)
  list_plots[[i]] = plot_grid(slog+theme(legend.position = "none"),
                              bplot+theme(legend.position = "none"),
                              ncol = 1, rel_heights = c(0.75,1), align= "v",axis = "lr")
  
}


scale = max(df_plot$group)
svg("HistogramFrustration_RAS_Rojasetal.svg", height = 4*scale, width = 15)
plot_grid(
  plot_grid(legend_logo, legend_bplot, ncol = 2) ,
  plot_grid(plotlist = list_plots,
            ncol = 1, align= "v",axis = "lr"),
  rel_heights = c(0.1,1.5),
  ncol = 1)
dev.off()

