#! /usr/bin/Rscript
################################################################

## LOAD PACKAGES
warn.conflicts = FALSE
options(echo = FALSE, verbose = F,warn = -1) 

requiredPackages = c('ggplot2', 'ggnewscale', 'cowplot',
                     'seqinr','stringr', 'data.table', "dplyr")#,'GetoptLong') 
suppressMessages(
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE))
      install.packages(p)
    library(p, character.only = TRUE)
  }
)

library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(data.table)
library(stringr)
library(dplyr)
############
### INPUT ##
############


#-------------------------------------------------------------------------
output = "./FrustraEvo_20221227610129084/pdb/FrustraEvo_20221227610129084"
seqrefs = "nsp3_protein_SARS-CoV-2"
cluster_names = "PLPro_cluster1"
msas = Sys.glob(file.path(output, "OutPutFiles/MSA_*.fasta"))
positions = Sys.glob(file.path(output, "AuxFiles/Positions"))
##########
### CODE##
##########

# extract and align single residue values

list_plots = list()
list_dfs = list()

for(m in 1:length(msas)){
  msa= read.fasta(msas[m],
                  forceDNAtolower = F, seqtype = "AA")
  
  positions_msa = readLines(positions[m])
  positions_msa = str_remove(positions_msa,'>')

  list_df = list()
  # get alignment and sequence indexes
  for (i in 1:length(msa)){
     if (str_detect(names(msa)[i], seqrefs[m])){
      IsRef = "yes"
    } else {
      IsRef = "no"
    }
    
    refid = names(msa)[i]
    
    n = which(positions_msa==refid)
    nref = which(positions_msa==seqrefs[m])
    positions_seq = data.frame ( id = positions_msa[c(n)],
                         Res = as.numeric(unlist(str_split(positions_msa[c(n+1)], ' '))),
                         Res_ref = as.numeric(unlist(str_split(positions_msa[c(nref+1)], ' '))))

    frust_sr = Sys.glob(file.path(output,'Data', paste(refid,"*done/*/*.pdb_singleresidue", sep ="")))
    cl = cluster_names[m]
    dt = fread(frust_sr)
    dt$Cluster= cl
    df2 = suppressMessages(full_join(dt, positions_seq))
    df2$Alignment_position = df2$Res_ref
    list_df[[i]] = df2
  }
  
  mapping_file = rbindlist(list_df)
  mapping_file = subset(mapping_file, !is.na(id))
  mapping_file$ id = factor(mapping_file$id,
                            levels = c(sort(unique(mapping_file$id)[-which(unique(mapping_file$id)==seqrefs[m])], decreasing = T),
                                       seqrefs[m]))
  
  
  mapping_file_subset = subset(mapping_file, Alignment_position %in% c(1:13,215:240))
  mapping_file_subset $g = ifelse(mapping_file_subset$Alignment_position %in% 1:13, '',
                                  ifelse(mapping_file_subset$Alignment_position%in% 215:220,"sel1",
                                         ifelse(mapping_file_subset$Alignment_position==221,"225",
                                                ifelse(mapping_file_subset$Alignment_position %in% 222:227,"sel2",
                                                       ifelse(mapping_file_subset$Alignment_position==228,"232","sel3")))))

  mapping_file_subset $g = factor( mapping_file_subset $g, levels= c("","sel1","225", "sel2","232", "sel3" ))
  ## PLOT
  g1 = ggplot(subset(mapping_file_subset, g == ""),
              aes(x=Alignment_position, y =id, fill = FrstIndex))+
    geom_tile( linejoin = "round",width = 0.96, height = 0.96, color = 'grey80')+
    geom_text(aes(label  = AA))+
    scale_fill_gradientn(name="Single Residue FI",
                         values = c(1, .79, .5, 0.1, 0), 
                         colours=c("green", "grey", "grey", "grey","red"),
                         limits=c(-1, 1),
                         oob = scales::squish,
                         na.value="white")+
    theme_minimal()+
    facet_grid(~g,drop = TRUE, scales="free_x")+
    ylab(NULL)+
    xlab(NULL)+
    theme(legend.position = "top")+
    coord_cartesian(ylim = c(0.25, length(unique(mapping_file$id)) + 0.5), clip="off") +
    scale_x_continuous( expand = c(0.001, 0.01),breaks = seq(0, 300, by = 10)) +
    scale_y_discrete(expand = c(0, 0)) 
  
 ellipsis =  ggplot(mapping_file_subset,
         aes(x=Alignment_position, y =id, fill = FrstIndex))+
    theme_classic()+
    annotate("text",label = '...', size = 16,y = 16, x = 3)+
   ylim(1,31)+
   xlim(1,5)#+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
 
    
    pos.labs <- c("", "225", "","232","")
    names(pos.labs) <- c("sel1","225","sel2", "232", "sel3")  
    
 g2 = ggplot(subset(mapping_file_subset, g != ""),
             aes(x=Alignment_position, y =id, fill = FrstIndex))+
   geom_tile( linejoin = "round",width = 0.96, height = 0.96, color = 'grey80')+
   geom_text(aes(label  = AA))+
   scale_fill_gradientn(name="Single Residue FI",
                        values = c(1, .79, .5, 0.1, 0), 
                        colours=c("green", "grey", "grey", "grey","red"),
                        limits=c(-1, 1),
                        oob = scales::squish,
                        na.value="white")+
   theme_minimal()+
   facet_grid(~g,drop = TRUE, scales="free_x", space="free_x",labeller = labeller(g=pos.labs))+
   ylab(NULL)+
   xlab(NULL)+
   theme(legend.position = "top")+
   coord_cartesian(ylim = c(0.25, length(unique(mapping_file$id)) + 0.5), clip="off") +
   scale_x_continuous( expand = c(0.001, 0.01),breaks = seq(0, 300, by = 10)) +
   scale_y_discrete(expand = c(0, 0)) 
 
 
  if(m == 1){
    legend= get_legend(g1)
    list_plots[[m]] = plot_grid(NULL,legend,
                                ncol =1,
                                rel_heights = c(0.5,1))
  }
  g1 = g1 +
    theme(legend.position = "none")
  
  g2 = g2 +
    theme(legend.position = "none",
          axis.text.y = element_blank())
  
  
  
  
  list_plots[[m+1]] = plot_grid(NULL, 
                                plot_grid(g1,
                                          ellipsis+theme(axis.title = element_blank(),
                                                         axis.text = element_blank(),
                                                         axis.ticks = element_blank(),
                                                         axis.line = element_blank()),
                                          g2,
                                          ncol = 3, rel_widths = c(0.9,0.20,1.4),
                                          align = "h"),
                                ncol =1,
                                rel_heights = c(0.05,1))
  
  list_dfs [[m]] = mapping_file
}

#### SAVE PLOT
m3 = rbindlist(list_dfs)
ns = unique(m3[, c("id","Cluster")])
h= c()
for (c in unique(ns$Cluster)){
  max_nseqs = max(table(ns$Cluster))
  sel = subset(ns, Cluster ==c)
  h = c(h, nrow(sel) / max_nseqs)
}


p = plot_grid(plotlist = list_plots,
              ncol =1,
              rel_heights = c(0.09, h),
              labels = c("", cluster_names))

library(grid)
library(gridExtra)
x.grob <- textGrob("Alignment position", 
                   gp=gpar( fontsize=12))

#add to plot

#p =grid.arrange(arrangeGrob(p, bottom = x.grob))

sc = 11
w = 108#unique(max(mapping_file$Alignment_position, na.rm=T))
#width = w/sc, height = 0.9*sc
svg('figureS10.svg',pointsize = "8",width = w/sc, height = 0.55*sc*length(cluster_names))
grid.arrange(arrangeGrob(p, bottom = x.grob))
invisible(dev.off())

