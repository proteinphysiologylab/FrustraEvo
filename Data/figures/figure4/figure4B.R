library(network)
library(sna)
library(ggplot2)
library(GGally)
library(data.table)
library(seqinr)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(stringr)
library(dplyr)
##################################
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
refpos$AA_res = paste(refpos$aa, refpos$RAS_ref_position, sep="")
mapping_file = merge(mapping_file, refpos[,c(2:3,6)], by="Alignment_position", allow.cartesian = T)

refs = data.frame(id = seqrefs, id2 = c( "ARF","RAS","RHO",  "RAB" ) )

i = c(10, 33, 53,65,68,75,81,83,115,145)
j = c(20, 37, 61,65,68,75,81, 83,121, 148)

v = c()
for (k in 1:length(i)){
  v = c(v,i[k]:j[k])
}

df = data.frame(res = v,
                domain = c(rep("G1", 11 ),
                           rep("G2",5),
                           rep("G3", 9),
                           rep('-',5),
                           rep("G4",7),
                           rep("G5",4)),
                sdp = c("no","yes",
                        rep ("no",2),"yes",
                        rep ("no",5),"yes",
                        "no","yes",
                        rep ("no",2),"yes",
                        rep ("no",3),"yes",
                        "no","yes","yes",
                        rep ("no",2),rep ("yes",5),
                        rep ("no",6),"yes",
                        rep ("no",4)))

mapping_file_motifs = merge(mapping_file, df, by.x ="RAS_ref_position", by.y= "res")
mapping_file_motifs_refs = subset(mapping_file_motifs, id %in% seqrefs)

####################
### only motifs and sdps
#####################
list_plots = list()
list_df = list()

for (prot in 1:length(seqrefs)){
  
  df3 = subset(mapping_file, id == seqrefs[prot])
  df3$AA_res = paste(df3$aa, df3$RAS_ref_position, sep ="")
  colnames(df3)[6] = "res"

  df4 = unique(merge(df3, df, by = "res"))

  ic_conf = fread(paste("./",
                        refs[prot,]$id2,
                        "/CMaps/IC_Mut", sep = ""),
                  header = T)
  colnames(ic_conf)[1:2] = c("Res1", "Res2")
  ic = subset(ic_conf,
               Res1 %in% as.character(df4$Protein_position) & 
                 Res2 %in% as.character(df4$Protein_position) & 
                ICtotal > 0.5 & EstadoConservado != "NEU" & FreqConts > 0.5)
  ic$Res1 = as.numeric(ic$Res1)
  ic$Res2 = as.numeric(ic$Res2)
  
  ic$color = ifelse(ic$EstadoConservado == "MIN", "green", ifelse(ic$EstadoConservado == "MAX", "red", "grey"))
  
  contacts = left_join(left_join(ic, df4, by = c("Res1"= "Protein_position") ),
            df4 , by = c("Res2"= "Protein_position"),suffix=c("1", "2") )
  contacts $ group =  refs[prot,]$id2
  
  list_df[[refs[prot,]$id2]] = contacts
  
  t = network(ic[,1:2])
  t %v% "col" = as.character(inner_join(data.frame(Protein_position= network.vertex.names(t)), df4)$domain)
  t %v% "names" = inner_join(data.frame(Protein_position= network.vertex.names(t)), df4)$AA_res
  t %v% "sdps" = as.character(inner_join(data.frame(Protein_position= network.vertex.names(t)), df4)$sdp)
  plot = ggnet2(t,
                edge.color = ic$color,
                label = "names",
                color = "col", 
                shape="sdps",
                size = 15, 
                label.size = 11,
                edge.size = 1.25,
                palette = c("G1" = "#e7d626",
                            "G2" = "#ffa317",
                            "G3" = "#fbb6dd",
                            "G4" = "#97f7fc",
                            "G5" = "#ff67ff",
                            "-" = "grey85"
                            
                ), 
                legend.size = 22,
                legend.position = "bottom", shape.legend = "SDP", color.legend = "Motif")+
    ggtitle(refs[prot,]$id2)+
    theme(title = element_text(size = 50, face = "bold"),
          legend.text = element_text(size=22 ))
  
  list_plots[[refs[prot,]$id2]] = plot
}

legend = get_legend(list_plots[[1]])

svg('figure4B.svg', width =42, height = 15)
plot_grid(plot_grid(list_plots[[1]]+theme(legend.position = 'none'),NULL,
                    list_plots[[2]]+theme(legend.position = 'none'),NULL,
                    list_plots[[3]]+theme(legend.position = 'none'),NULL,
                    list_plots[[4]]+theme(legend.position = 'none'),
                    ncol = 7  ,
                    rel_widths = c(1.1,0.1,1.1,0.1,1.2,0.1,1)),
          legend,
          ncol = 1,
          rel_heights = c(1,0.2))
dev.off()



