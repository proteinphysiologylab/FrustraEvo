#! /usr/bin/Rscript
################################################################

## LOAD PACKAGES
warn.conflicts = FALSE
options(echo = FALSE, verbose = F,warn = -1) 

#also required Biostrings!!!
requiredPackages = c('ggplot2', 'ggnewscale', 'cowplot',
                     'seqinr','stringr', 'data.table', "dplyr", "tidyr","ggdendro","ggtext",'GetoptLong') 
suppressMessages(
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE))
      install.packages(p)
    library(p, character.only = TRUE)
  }
)

############
### INPUT ##
############

GetoptLong(
  #"msas=s@", "List of MSA fasta files",
  "output=s@", "List of FrustraEvo output folders",
  "seqrefs=s@", "List of sequence reference IDs.",
  "cluster_names=s@", "Cluster names of each MSA file."
)

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
                                 Res_seq = as.numeric(unlist(str_split(positions_msa[c(n+1)], ' '))),
                                 Res_ref = as.numeric(unlist(str_split(positions_msa[c(nref+1)], ' '))))
    
    v = as.numeric(unlist(str_split(positions_msa[c(n+1)], ' ')))
    v = v[-length(v)]
    
    frust_sr = Sys.glob(file.path(output,'Data', paste(refid,".done/*/*.pdb_singleresidue", sep ="")))
    cl = cluster_names[m]
    dt = fread(frust_sr)
    dt$Cluster= cl
    
    
    if (length(dt$Res_pdb) != length(v)){
      colnames(dt)[colnames(dt)=="Res"] = "Res_pdb"
      library("Biostrings")
      seq=c2s(getSequence(msa[i])[[1]][getSequence(msa[i])[[1]] != "-"])
      pdb_seq =c2s(dt$AA)
      PWA  = pairwiseAlignment(seq, seq)
      
      posss= data.frame (AA_seq = s2c(as.character(pattern(PWA))),
                         AA_pdb = s2c(as.character(subject(PWA))))
      
      start = max(PWA@subject@range@start,PWA@pattern@range@start)
      Res_seq = start
      Res_seq_v=c()
      
      Res_pdb = start
      Res_pdb_v=c()
      for (j in 1:nrow(posss)){
        
        if (posss$AA_seq[j] !="-"){
          
          Res_seq_v = c(Res_seq_v, Res_seq )
          Res_seq = Res_seq+1
        } else{
          Res_seq_v = c(Res_seq_v, NA )
        }
        
        if (posss$AA_pdb[j] !="-"){
          
          Res_pdb_v = c(Res_pdb_v, Res_pdb )
          Res_pdb = Res_pdb+1
        } else{
          Res_pdb_v = c(Res_pdb_v, NA )
        }
        
      }
      posss$Res_seq = Res_seq_v
      posss$Res_pdb = Res_pdb_v
      
      poss_filtered = suppressMessages(left_join(positions_seq,posss))
      df2 = left_join( poss_filtered, dt, by = "Res_pdb")
      df2$FrstIndex = ifelse(df2$AA_seq=="-", NA, df2$FrstIndex)
    }else{
      colnames(dt)[colnames(dt)=="Res"] = "Res_seq"
      df2 = left_join(positions_seq,dt, by ="Res_seq")
    }
    
    
    df2$Alignment_position = df2$Res_ref
    list_df[[i]] = df2
  }
  
  mapping_file = rbindlist(list_df, fill = T)
  
  # preprocessing for dendrogram    
  toreshape = mapping_file[,c("Alignment_position","FrstIndex","id")]
  tore <- toreshape# %>% distinct
#  rs = pivot_wider(toreshape, names_from = Alignment_position, values_from = FrstIndex)
  rs = spread(tore, Alignment_position, FrstIndex)
  rs$`<NA>` = NULL
  rsm = as.matrix(rs[,-c("id")])
  rownames(rsm) = rs$id
  # calculate clusters' dendrogram based on euclidian distance
  dend = as.dendrogram(hclust(dist(as.matrix(rsm),method = "manhattan")))
  dend_data <- dendro_data(dend)
  
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  
  pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, axpos = as.character(label), height = 1))  
  
  axis_limits <- with(
    pos_table, 
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
  ) + 
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  mapping_file$id = factor(mapping_file$id,
                           levels = dend_data$labels$label)
  #levels = c(sort(unique(mapping_file$id)[-which(unique(mapping_file$id)==seqrefs[m])], decreasing = T),
  #           seqrefs[m]))
  
  bold.labels <- ifelse(levels(mapping_file$id) %in% seqrefs, yes = paste("**",seqrefs,"**", sep =""),
                        no = levels(mapping_file$id))
  
  plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_x_reverse(expand = c(0, 0.5)) + 
    scale_y_continuous(breaks = pos_table$y_center, 
                       labels = bold.labels, 
                       limits = axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_markdown())
  
  
  
  ## PLOT
  g1 = ggplot(mapping_file,
              aes(x=Alignment_position, y =id, fill = FrstIndex))+
    geom_tile( linejoin = "round",width = 0.96, height = 0.96, color = 'grey80')+
    geom_text(aes(label  = AA_seq))+
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
    scale_y_discrete(expand = c(0, 0),labels = bold.labels) 
  
  
  if(m == 1){
    legend= get_legend(g1)
    list_plots[[m]] = plot_grid(NULL,NULL,NULL,legend,
                                ncol =2,
                                rel_heights = c(0.5,1))
  }
  g1 = g1 +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = unit(c(1, 0.2, 0.2, 0), "cm"))
  
  list_plots[[m+1]] = plot_grid(NULL, NULL, plt_dendr,g1,
                                ncol =2,
                                rel_heights = c(0.05,1),
                                rel_widths = c(0.20,1),
                                align = "h")
  
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
              rel_heights = c(1, max_nseqs),
              labels = c("", cluster_names))
sc = 7
a=length(positions_msa)/42
w = unique(max(mapping_file$Alignment_position, na.rm=T))
#width = w/sc, height = 0.9*sc
svg('MSFA.svg',pointsize = "100",idth = 100, height = 30)
p
invisible(dev.off())
library(magick)
testimage <- image_read_svg('MSFA.svg', width = 2500)
image_write(testimage, path = "MSFA.png", format = "png")
