library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(tidyr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
##################################
icfrst = Sys.glob("./IC_seqfrust.csv")
icfrst_f = list()
for (i in icfrst){
  cl = "RAS"
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

### CORRelation
fitness= read_xlsx('../media-4.xlsx',sheet = 2)
colnames(fitness)

unique(fitness$assay)
unique(subset(fitness, WT ==TRUE)$aa_seq)

fitness$Protein_position = mapply(function(x, y) which(x != y)[1], 
       strsplit(fitness$aa_seq[1], ""), strsplit(fitness$aa_seq, ""))

head(fitness$Protein_position)

fitness$Protein_position = fitness$Protein_position + 1

ictot_icseq_fitness = merge(ictot_icseq, fitness, by ="Protein_position", allow.cartesian = T)


ictot_icseq_fitness$Assay = str_split_fixed(ictot_icseq_fitness$assay, ' ', 2)[,1]

df_plot = ictot_icseq_fitness%>% group_by(Protein_position,Assay,  ICTot, ICseq, FrustEstado)%>%
   summarise(mean_fitness = mean(fitness), 
             mean_sigma = mean(sigma))
 

s1= ggscatter(subset(df_plot, Assay == "AbundancePCA"),
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              add = "reg.line",
              palette= c("red","green3", "grey"),
              ylab = "FrustIC",
              xlab="Abundance fitness",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5), xlim = c(-1,0.25))+
  stat_cor(aes(color = FrustEstado),
           label.x = -0.35,label.y = c(1.25,1.15,1.05),
           size = 3.5)
s1

saveRDS(s1, './KRAS_corr_icfrust_vs_abundancefitness.RDS')


s1= ggscatter(subset(df_plot, Assay != "AbundancePCA"),
              x = "mean_fitness",
              y = "ICTot",
              color = "FrustEstado",
              # cor.coef = T, 
              add = "reg.line",
              palette= c("red","green3", "grey"),
              ylab = "FrustIC",
              xlab="Binding fitness",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
  stat_cor(aes(color = FrustEstado),
           label.x = -1.3, 
           size = 3.5)
s1

saveRDS(s1, './KRAS_corr_icfrust_vs_bindingfitness.RDS')



s2= ggscatter(subset(df_plot, Assay == "AbundancePCA"),
              x = "mean_fitness",
              y = "ICseq",
              color = "FrustEstado",
               cor.coef = T, 
              cor.coef.size = 3.5,
              cor.coef.coord = c(-0.35,3.75),
              #add = "reg.line",
              palette= c("red","green3", "grey"),
              ylab = "SeqIC",
              xlab="Abundance fitness",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
  geom_smooth(method="lm", linewidth = 0.5)
s2

saveRDS(s2, './KRAS_corr_icseq_vs_abundancefitness.RDS')


s2= ggscatter(subset(df_plot, Assay != "AbundancePCA"),
              x = "mean_fitness",
              y = "ICseq",
              color = "FrustEstado",
               cor.coef = T, 
              cor.coef.size = 3.5,
              palette= c("red","green3", "grey"),
              ylab = "SeqIC",
              xlab="Binding fitness",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
  geom_smooth(method="lm", linewidth = 0.5)
s2

saveRDS(s2, './KRAS_corr_icseq_vs_bindingfitness.RDS')

s3= ggscatter(ictot_icseq,
              x = "ICseq",
              y = "ICTot",
              color = "FrustEstado",
              cor.coef.size = 3.5,
              # cor.coef = T, 
              add = "reg.line",
              palette= c("red","green3", "grey"),
              ylab = "FrustIC",
              xlab="SeqIC",
              font.label = c(14, "plain"),
              size = 1,
              add.params = list(size = 0.5))+
  stat_cor(aes(color = FrustEstado), label.x = 0, size = 3.5)
s3
saveRDS(s3, './KRAS_corr_icseq_vs_icfrust.RDS')
