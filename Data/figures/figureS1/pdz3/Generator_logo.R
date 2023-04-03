# load required packages
suppressMessages(require(ggplot2))
suppressMessages(require(ggseqlogo))
suppressMessages(library(cowplot))
suppressMessages(library(seqinr))
suppressMessages(library(data.table))
library(ggnewscale)
library(dplyr)

# load msa to generate logo
msa = read.fasta(file = "./SalidaAlignSE.fasta", forceDNAtolower = F, seqtype = "AA")

# read frust data
FileTabla2=read.table(file="./CharactPosDataN", header=T, stringsAsFactors=F)

h_minimos=as.numeric(FileTabla2[,8])
h_neutros=as.numeric(FileTabla2[,9])
h_maximos=as.numeric(FileTabla2[,10])

N=length(h_minimos)

#generate data for bars
df_plot = data.table(N=rep(1:N,3),
                       h= c(h_minimos, h_maximos,h_neutros),
                       state= rep(c('Minimum', "Maximum", "Neutral"), each=N),
                       group=1)
start = 15
end = 98


load("./JD_PDZ_NM2_stabilityPCA_dimsum128_filtered_fitness_replicates.RData")

fitness = singles
fitness$Protein_position = fitness$Pos + start -1
fitness2 = fitness%>%
  group_by(Protein_position)%>%
  summarise(median_fitness = median(fitness))

df_plot = inner_join(df_plot, fitness2, by = c("N"="Protein_position"))


# generate logo
slog = ggseqlogo(unlist(getSequence(getFrag(msa, start,end), as.string = T)),
                 ncol = 100,
                 font= "akrobat_bold")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y =  element_text(size = 16))+
  scale_x_continuous(breaks = c(1,seq(5,2998, by = 5)), expand = c(0, 0))

LOGO_DATA = ggplot_build(slog)$data[[1]]
LOGO_DATA$roundx = round(LOGO_DATA$x)
max = LOGO_DATA %>% group_by(roundx) %>% summarise(max = max(y))
avg_entropy = - (mean(max$max) - log2(20))


# Plot
slog = slog + 
  ylab("SeqIC")

max_ic_frust = df_plot %>% group_by(N) %>% summarise(max = max (h))
avg_entropy_frust = - (mean(max_ic_frust $max) - log2(3))

# generate barplot
bplot = ggplot(data = subset(df_plot, group == 1 ), aes(x=N, y =h, fill = state))+
    geom_bar(stat= "identity",position="stack", color = "black")+
    scale_fill_manual(values= c("red4", "darkgreen", "grey"))+
    scale_x_continuous(breaks = c(1,seq(5,1998, by = 5)), expand = c(0, 0))+
   theme_logo()+
    xlab("")+
    ylab("FrustIC")+
    theme(axis.text = element_text(size = 12),
          axis.ticks = element_line(),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.y =  element_text(size = 16))


bplot
  
  legend_logo = get_legend(slog)
  legend_bplot = get_legend(bplot)
   plot_grid(slog+theme(legend.position = "none"),
                              bplot+theme(legend.position = "none"),
                              ncol = 1, rel_heights = c(0.75,1), align= "v",axis = "lr")
  


# save figure
scale = max(df_plot$group)
svg("HistogramFrustration_PDZ.svg", height = 4*scale, width = 15)
plot_grid(slog+theme(legend.position = "none"),
          bplot+theme(legend.position = "none"),
          ncol = 1, rel_heights = c(0.75,1), align= "v",axis = "lr")
dev.off()
