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
# LOGO PSDZ

msa = read.fasta('./filtered_renamed_PSD95.fasta',
                 forceDNAtolower = F, seqtype = "AA")


seqrefs= c("1BE9_A")


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
files2$Protein_position =  files2$Res
files2
ictot_icseq = files2
ictot_icseq$ICTot = as.numeric(ictot_icseq$ICTot)
ictot_icseq$ICseq = as.numeric(ictot_icseq$ICseq)


i = 15
j = 98


##### BINDING #######
load("./JD_PDZ_NM2_bindingPCA_dimsum128_filtered_fitness_replicates.RData")

fitness = singles
fitness$Protein_position = fitness$Pos + i -1

ictot_icseq_fitness = merge(ictot_icseq, fitness, by ="Protein_position", allow.cartesian = T)

df_plot = ictot_icseq_fitness%>% group_by(Pos,Cluster,  ICTot, ICseq, FrustEstado, WT_AA)%>%
   summarise(mean_fitness = mean(fitness), 
             mean_sigma = mean(sigma))
 
df_plot$label = paste(df_plot$WT_AA, df_plot$Pos, sep="")



s1= ggscatter(df_plot,
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              size = 1,
              # cor.coef = T, 
              add = "reg.line",
              add.params = list(size=0.5),
              palette= c("green3", "grey"),
              ylab = "FrustIC",
              xlab="Binding fitness",
              font.label = c(14, "plain"))+
  stat_cor(aes(color = FrustEstado), 
           label.x = -0.5,
           size =3.5)


saveRDS(s1, './PSD95_corr_icfrust_vs_bindingfitness.RDS')

s2= ggscatter(df_plot,
              x = "mean_fitness",
              y = "ICseq",
              color = "FrustEstado",
              cor.coef = T,
              cor.coef.coord = c(-0.6,4) ,
              cor.coef.size = 3.5,
              size = 1,
              #add = "reg.line",
              palette= c("green3", "grey"),
              ylab = "SeqIC",
              xlab="Binding fitness",
              font.label = c(14, "plain"))+
  geom_smooth(method="lm", linewidth = 0.5)



s2
saveRDS(s2, './PSD95_corr_icseq_vs_bindingfitness.RDS')


s3= ggscatter(df_plot,
              x = "ICseq",
              y = "ICTot",
              color = "FrustEstado",
              # cor.coef = T, 
              add = "reg.line",
              add.params = list(size = 0.5),
              palette= c("green3", "grey"),
              ylab = "FrustIC",
              xlab="SeqIC",
              size = 1,
              font.label = c(14, "plain"))+
  stat_cor(aes(color = FrustEstado), label.x = 1, size = 3.5)
s3
saveRDS(s3, './PSD95_corr_icfrust_vs_icseq.RDS')


### STABILITY #######
load("./JD_PDZ_NM2_stabilityPCA_dimsum128_filtered_fitness_replicates.RData")

fitness = singles
fitness$Protein_position = fitness$Pos + i -1

ictot_icseq_fitness = merge(ictot_icseq, fitness,allow.cartesian = T)

df_plot = ictot_icseq_fitness%>% group_by(Pos,Cluster, ICTot, ICseq, FrustEstado)%>%
   summarise(mean_fitness = mean(fitness), 
             mean_sigma = mean(sigma))

s1= ggscatter(df_plot,
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              # cor.coef = T, 
              add = "reg.line",
              add.params = list(size = 0.5),
              palette= c("green3", "grey"),
              ylab = "FrustIC",
              xlab = "Abundance fitness",
              size = 1,
              font.label = c(14, "plain"))+
  stat_cor(aes(color = FrustEstado),
           label.x = -0.45,
           size = 3.5)
s1
saveRDS(s1, './PSD95_corr_icfrust_vs_abundancefitness.RDS')

s2=  ggscatter(df_plot,
               x = "mean_fitness",
               y = "ICseq",
               color = "FrustEstado",
               cor.coef = T,cor.coef.coord = c(-0.4,3.8) ,
               cor.coef.size = 3,
               palette= c("green3", "grey"),
               ylab = "SeqIC",
               xlab = "Abundance fitness",
               size = 1,
               font.label = c(14, "plain"))+
  geom_smooth(method="lm",linewidth = 0.5)

s2
saveRDS(s2, './PSD95_corr_icseq_vs_abundancefitness.RDS')

