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

output = "./FrustraEvo_Alfas"
seqrefs = "2dn1-A"
cluster_names = ""
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
  
  
  mapping_file_subset = subset(mapping_file, Alignment_position %in% c(1:55))
  ## PLOT
  g1 = ggplot(subset(mapping_file_subset),
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
  
  list_plots[[m+1]] = plot_grid(NULL, 
                                g1,
                                ncol =1,
                                rel_heights = c(0.05,1))
  
  list_dfs [[m]] = mapping_file_subset
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


sc = 11
w = 10
svg('MSFA.svg',pointsize = "8",width = w/sc, height = 0.55*sc*length(cluster_names))
p
invisible(dev.off())

