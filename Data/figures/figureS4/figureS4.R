library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
library(grid)
library(gridExtra)
################################################
## DECOMPRESS FILES IN THE DIRECTORY FIRST!!!!
################################################
msa = read.fasta('./ras_human.fasta',
                 forceDNAtolower = F, seqtype = "AA")

names(msa) = str_split_fixed(names(msa),pattern = '\\|',n = 3)[,2]
names(msa) = str_remove(names(msa), '-1')

seqrefs= c("P84077" ,"P01112", "P61586", "P62820")

list_df = list()

for (i in 1:length(msa)){
  s = msa[[i]]
  index_a = 0
  index_s = 0
  aa = c()
  seqpos = c()
  alipos = c()
  for(j in 1:length(s)){
    if (s[j] !='-'){
      index_a = index_a + 1
      index_s = index_s + 1
      aa = c(aa, s[j])
      seqpos = c(seqpos, index_s)
      alipos = c(alipos, index_a)
    } else {
      index_a = index_a + 1
    }
  }
  if (names(msa)[i] %in% seqrefs){
    IsRef = "yes"
  } else {
    IsRef = "no"
  }
  
  df = data.frame(aa, Protein_position = seqpos, Alignment_position = alipos, id= names(msa)[i], IsRef )
  list_df[[i]] = df
}

mapping_file = rbindlist(list_df)
mapping_file

refpos = subset(mapping_file, id == "P01112")
colnames(refpos)[2] = "RAS_ref_position"
mapping_file = merge(mapping_file, refpos[,2:3], by="Alignment_position", allow.cartesian = T)

cl1 = Sys.glob("./*/OutPutFiles*/Frustration/*done/FrustrationData/*.pdb_singleresidue")
library(stringr)
cl1_f = list()
for (i in cl1){
  id = str_remove(basename(dirname(dirname(i))), '\\.done')
  cl = basename(dirname(dirname(dirname(dirname(dirname(i))))))
  dt = fread(i)
  dt$id = id
  dt$Cluster= cl
  cl1_f [[i]] = dt
}
files = rbindlist(cl1_f)
colnames(files)[1] = "Protein_position"

m2 = merge(files, mapping_file , by = c("Protein_position", "id"))

m2 = m2 %>% group_by (Cluster, Alignment_position)%>%
  mutate(median_FrstIndex = median(FrstIndex))

mref = subset(m2, IsRef == "yes")

icfrst = Sys.glob("./*/OutPutFiles*/Equivalences/IC_seqfrust.csv")
icfrst_f = list()
for (i in icfrst){
  cl = str_remove(basename(dirname(dirname(dirname(i)))), 'cluster')
  dt = fread(i, dec = ",")
  dt$Cluster= cl
  dt$ICseq = as.numeric(dt$ICseq)
  icfrst_f [[i]] = dt
}
files2 = rbindlist(icfrst_f)
files2
files2$Protein_position =  files2$Res

ictot_icseq = merge(mref, files2)


i = c(10, 33, 53,65,68,75,81,83,115,145)
j = c(20, 37, 61,65,68,75,81, 83,121, 148)

ictot_icseq$Cluster = factor(ictot_icseq$Cluster, levels = c("RAB","RHO","RAS","ARF"))
ictot_icseq$threshold = ifelse(ictot_icseq$ICTot > 0.5, 'ICTot > 0.5', 'ICTot <= 0.5')
ictot_icseq$median_FrstIndex2  = ictot_icseq$median_FrstIndex
ictot_icseq[ictot_icseq$ICTot <= 0.5,]$median_FrstIndex2 = NA
ictot_icseq$label = ifelse(ictot_icseq$RAS_ref_position %in% c(11,14,20,34,37,56,58,59,65,68,69,75,81,83,121), '*', NA)

list_plots = list()
list_plots_ann = list()
for (k in 1:length(i)){
  s = subset(ictot_icseq, RAS_ref_position %in% i[k]:j[k] )

  g1 = ggplot(s,
              aes(x= as.character(RAS_ref_position), y = Cluster, fill = median_FrstIndex2))+
    geom_tile( linejoin = "round",width = 0.96, height = 0.96, color = 'grey80')+
    geom_text(aes(label  = AA, size = ICseq), fontface = "bold")+
    scale_fill_gradientn(values = c(1, .79, .5, 0.1, 0), 
                         colours=c("green", "grey", "grey", "grey","red"),
                         limits=c(-1, 1),
                         oob = scales::squish,
                         na.value = "white",
                         name= "median SRFI")+
    scale_color_gradient2(mid = "grey", high = "black",midpoint=0.5, limits=c(0, log2(20)),
                          oob = scales::squish,
                          name= "median SRFI")+
    theme_minimal()+
    ylab(NULL)+
    xlab(NULL)+
    theme(legend.position = "none",
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.text.x = element_text(size=12))+
    coord_cartesian(ylim = c(0.5, 4.5), clip="off") +
    scale_y_discrete(expand = c(0, 0)) 
  
  if (k != 1){
    g1 = g1+theme(axis.text.y = element_blank())
  }
  ann <- ggplot(s, aes(x = as.character(RAS_ref_position), y = 0, label =label)) +
    geom_text(size = 10) +
    theme_void()
  
  list_plots[[k]] = g1
  list_plots_ann[[k]] = ann
}
  
# combine both plot using plot_grid()
combined_plot<-plot_grid(
  plot_grid(NULL, list_plots_ann[[1]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[1]]+
              theme(panel.border = element_rect(color = "#e7d626", fill=NA, linewidth = 1.5)),
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("               G1", "",""),  label_size = 20,align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[2]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[2]]+
              theme(panel.border = element_rect(color = "#ffa317", fill=NA, linewidth = 1.5)),
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("     G2", "",""), label_size = 20, align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[3]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[3]]+
              theme(panel.border = element_rect(color = "#fbb6dd", fill=NA, linewidth = 1.5)),
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("        G3", "",""),  label_size = 20,align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[4]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[4]],
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("", "",""), align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[5]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[5]],
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("", "",""), align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[6]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[6]],
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("", "",""), align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[7]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[7]],
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("", "",""), align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[8]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[8]],
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("", "",""), align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[9]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[9]]+
              theme(panel.border = element_rect(color = "#97f7fc", fill=NA, linewidth = 1.5)),
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("          G4", "",""), label_size = 20, align = "v", axis = "l"),
  plot_grid(NULL, list_plots_ann[[10]]+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            list_plots[[10]]+
              theme(panel.border = element_rect(color = "#ff67ff", fill=NA, linewidth = 1.5)),
            ncol = 1, rel_heights = c(0.2,0.1,1), labels = c("    G5", "","") ,label_size = 20, align = "v", axis = "l"),
  ncol=10,
  rel_widths = c(1,0.5,0.75,0.1,0.1,0.1,0.1,0.1,0.75,0.5))

combined_plot

legend2 <- get_legend(
  g1 +
    guides(
           size = guide_legend("SeqIC")) +
    theme(legend.position = "right")
)

# Combine combined plot and legend using plot_grid()
p = plot_grid(  combined_plot,legend2, ncol=2,rel_widths =  c(1, .1))


y.grob <- textGrob("", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13), rot=90)

x.grob <- textGrob("RAS ref position", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13))


svg('figureS4.svg', 
    height = 4,
    width = 16)
grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))
dev.off()

