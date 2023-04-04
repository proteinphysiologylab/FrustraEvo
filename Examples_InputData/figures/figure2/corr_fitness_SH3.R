library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggpubr)
##################################
msa = read.fasta('./filtered_renamed_GRB2.fasta',
                 forceDNAtolower = F, seqtype = "AA")

icfrst = Sys.glob("./IC_seqfrust.csv")
icfrst_f = list()
for (i in icfrst){
  cl = basename(dirname(dirname(i)))
  dt = fread(i, dec = ",")
  dt$Cluster= cl
  icfrst_f [[i]] = dt
}
files2 = rbindlist(icfrst_f)
files2
files2$Ref_position =  as.numeric(files2$Res)
files2$Ref_position = files2$Ref_position -2

ictot_icseq = files2
ictot_icseq$ICTot = as.numeric(ictot_icseq$ICTot)
ictot_icseq$ICseq = as.numeric(ictot_icseq$ICseq)

#### BINDING
load("./JD_GRB2_epPCA_bindingPCA_dimsum128_fitness_replicates.RData")

fitness = singles
fitness$Ref_position = fitness$Pos 

ictot_icseq_fitness = merge(ictot_icseq, fitness, by ="Ref_position", allow.cartesian = T)

df_plot = ictot_icseq_fitness%>% group_by(Pos, Cluster,ICTot, ICseq, FrustEstado, WT_AA)%>%
  summarise(mean_fitness = mean(fitness), 
            mean_sigma = mean(sigma))

df_plot$label = paste(df_plot$WT_AA, df_plot$Pos, sep="")

s1= ggscatter(df_plot,
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              # cor.coef = T, 
              add = "reg.line",
              palette= c("red", "green3", "grey"),
              ylab = "FrustIC",
              xlab="Binding fitness",
              size=1,
              add.params = list(size = 0.5),
              font.label = c(14, "plain"))+
  stat_cor(aes(color = FrustEstado),
           label.x = -0.45, size = 3.5)
s1
saveRDS(s1, './SH3_corr_icfrust_vs_bindingfitness.RDS')
s2=  ggscatter(df_plot,
               x = "mean_fitness",
               y = "ICseq",
               color = "FrustEstado",
               cor.coef = T,cor.coef.coord = c(-0.5,4) ,
               cor.coef.size = 3.5,
               #add = "reg.line",
               palette= c("red","green3", "grey"),
               ylab = "SeqIC",
               xlab="Binding fitness",
               font.label = c(14, "plain"),
               size = 1,
               add.params = list(size = 0.5))+
  geom_smooth(method="lm", linewidth=0.5)

s2
saveRDS(s2, './SH3_corr_icseq_vs_bindingfitness.RDS')

s3= ggscatter(df_plot,
              x = "ICseq",
              y = "ICTot",
              color = "FrustEstado",
              cor.coef = F, 
              add = "reg.line",
              palette= c("red", "green3", "grey"),
              ylab = "FrustIC",
              xlab = "SeqIC",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
 stat_cor(aes(color = FrustEstado),
          label.x = 0.5,
          size = 3.5)#+
  #geom_smooth(method="lm")

s3
saveRDS(s3,'./SH3_corr_icseq_vs_icfrust.RDS')



#### STABILITY
load("./JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData")

fitness = singles
fitness$Ref_position = fitness$Pos 

ictot_icseq_fitness = merge(ictot_icseq, fitness)

df_plot = ictot_icseq_fitness%>% group_by(Pos, ICTot, ICseq, FrustEstado,WT_AA)%>%
  summarise(mean_fitness = mean(fitness), 
            mean_sigma = mean(sigma))
df_plot$label = paste(df_plot$AA, df_plot$Pos, sep="")

s1= ggscatter(df_plot,
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              # cor.coef = T, 
              add = "reg.line",
              palette= c("red", "green3", "grey"),
              ylab = "FrustIC",
              xlab="Abundance fitness",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
  stat_cor(aes(color = FrustEstado), label.x = -0.5, size = 3.5)

s1
saveRDS(s1,'./SH3_corr_icfrust_vs_abundancefitness.RDS')


s2=  ggscatter(df_plot,
               x = "mean_fitness",
               y = "ICseq",
               color = "FrustEstado",
               cor.coef = T,cor.coef.coord = c(-0.5,4) ,
               cor.coef.size = 3.5,
               #add = "reg.line",
               palette= c('red', "green3", "grey"),
               ylab = "SeqIC",
               xlab="Abundance fitness",
               font.label = c(14, "plain"),
               size = 1,
               add.params = list(size = 0.5))+
  geom_smooth(method="lm", linewidth = 0.5)
s2
saveRDS(s2,'./SH3_corr_icseq_vs_abundancefitness.RDS')
