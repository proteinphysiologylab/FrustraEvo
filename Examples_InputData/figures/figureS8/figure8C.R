library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
library(ggpubr)
library(viridis)
##################################

wd = "../../figure5/SARS2_FrustraEvo_input_and_output"

prot_df = data.frame(
  protein = c(
    'nsp1',
    'nsp2',
    'nsp3_Ubl1',
    'nsp3_MacroX',
    'nsp3_MacroII',
    'nsp3_SUD_Nterm',
    'nsp3_SUD_Cterm',
    'nsp3_PLPro',
    'nsp3_NAB',
    'nsp3_Y3',
    'nsp5',
    'nsp7',
    'nsp8',
    'nsp9',
    'nsp10',
    'nsp12',
    'nsp13',
    'nsp14',
    'nsp15',
    'nsp16',
    'E',
    'N_Cterm',
    'N_Nterm',
    'S_RBD',
    "S_protein",
    'orf3a',
    'orf7a',
    'orf8',
    'orf9b'
  ),
  old_name = c(
    'nsp1',
    'nsp2',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp3',
    'nsp5',
    'nsp7',
    'nsp8',
    'nsp9',
    'nsp10',
    'nsp12',
    'nsp13',
    'nsp14',
    'nsp15',
    'nsp16',
    'E',
    'N',
    'N',
    'S',
    "S",
    'orf3a',
    'orf7a',
    'orf8',
    'orf9b'
  ),
  
  prot = c(
    'NSP1_protein',
    'NSP2_protein',
    'NSP3_cl13772_DUF3655',
    'NSP3_cd21557_Macro_X_Nsp3-like',
    'NSP3_cl00019_Macro_SF',
    'NSP3_cl13138_SUD-M',
    'NSP3_cd21525_SUD_C_SARS-CoV_Nsp3',
    'NSP3_cd21732_betaCoV_PLPro',
    'NSP3_cd21822_SARS-CoV-like_Nsp3_NAB',
    'NSP3_cd21717_TM_Y_SARS-CoV-like_Nsp3_C',
    'NSP5_protein',
    'NSP7_protein',
    'NSP8_protein',
    'NSP9_protein',
    'NSP10_protein',
    'NSP12_protein',
    'NSP13_protein',
    'NSP14_protein',
    'NSP15_protein',
    'NSP16_protein',
    'E_protein',
    'N_C-terminaldomain',
    'N_N-terminaldomain',
    'S_RBD',
    "S_protein",
    'orf3a_protein',
    'orf7a_protein',
    'orf8_protein',
    'orf9b_protein'
  )
)

select= c(   'NSP2_protein',
             'NSP3_cl13772_DUF3655',
             'NSP3_cd21557_Macro_X_Nsp3-like',
             'NSP3_cl13138_SUD-M',
             'NSP3_cd21525_SUD_C_SARS-CoV_Nsp3',
             'NSP3_cd21732_betaCoV_PLPro',
             'NSP3_cd21822_SARS-CoV-like_Nsp3_NAB',
             'NSP3_cd21717_TM_Y_SARS-CoV-like_Nsp3_C',
             'NSP5_protein',
             'NSP7_protein',
             'NSP8_protein',
             'NSP9_protein',
             'NSP10_protein',
             'NSP12_protein',
             'NSP13_protein',
             'NSP14_protein',
             'NSP15_protein',
             'NSP16_protein',
             'E_protein',
             'N_C-terminaldomain',
             'N_N-terminaldomain',
             "S_protein")



nseqs_plddt = fread('../../figureS6/nseqs_plddt.txt')
select = subset(nseqs_plddt, nseqs>=10 & y >= 80)
colnames(select)[colnames(select)=="cluster"] = "Cluster"


taxid = fread("./taxid_sdpsclusters.txt")
subset_taxid = unique(taxid[,c("seqids","taxids","cl_phylo","subcl_phylo")])
colnames(subset_taxid)[1] = "Feature"
########
# SDPs vs Frust 
########
mapper_fs =  Sys.glob(file.path(wd, "*/*mapped_positions_unique.csv"))

# read 3dmapper results and merge them with frustration single res results
pident_results = list()
frust_results = list()
for(f in mapper_fs){
  pid =basename(dirname(f))
  m = fread(f)
  pids = str_split_fixed(unique(m$Feature), '\\.',2)[,1]
  blast_results = fread(Sys.glob(paste("./blast_raw_results/",
                              subset(prot_df, prot %in% pid)$old_name,
                              ".blast",
                        sep ="")))
  
  subset_blast= subset(blast_results,
                        V4 %in% pids)[,c("V1","V3","V4","V5","V13")]
  subset_blast$seqids = str_split_fixed(subset_blast$V3,"\\|", 3)[,2]
  
  colnames(subset_blast) = c("query", 'seqid',"accid", "taxids", "pident","Feature")
  
  subset_blast= (merge(subset_blast, unique(m[,c("Gene", "Feature","Cluster")])))
  colnames(subset_blast)[colnames(subset_blast)=="Gene"] = "prot"
  
  subset_blast = merge(subset_blast, prot_df, by = "prot")
  
  subset_blast$taxids = as.numeric(subset_blast$taxids)
  
  subset_blast = merge(subset_blast,subset_taxid, by = c("taxids", "Feature"))
  pident_results[[pid]] = subset_blast
}


pident_all = rbindlist(pident_results)
pident_all$n = 1
pident_all$cl_phylo[pident_all$cl_phylo=="Unclassified"]="Betacoronavirus"

pident_all_sum = unique(pident_all %>% group_by (prot, Cluster) %>%
  mutate(total = sum (n))%>%
  group_by (prot, Cluster, subcl_phylo, cl_phylo) %>%
  summarise(dens = sum (n)/total))

pident_all_sum$n = 1

subgenus_count = pident_all_sum%>% group_by(prot, Cluster, cl_phylo) %>% summarise(n_subgenus = sum(n), sd_dens = sd(dens))
subgenus_count$sd_dens = ifelse(is.na(subgenus_count$sd_dens), 0,subgenus_count$sd_dens)


res_df = fread('../../figureS7/frust_3dmapper_sdps_results.txt')
colnames(res_df)[colnames(res_df)=="Gene"] = "prot"
res_df_sum =subset(res_df, subcluster == "rest") %>% group_by(prot, Cluster) %>%
  summarise(mean_icseq= mean(ICseq), mean_icfrust = mean(ICTot))

res_df_sum_filtered = merge(res_df_sum, select)


res_df_sum_filtered = left_join(res_df_sum_filtered, subgenus_count, by =c("prot","Cluster"))


scatter = ggscatter(subset(res_df_sum_filtered, prot !="S_RBD" & sd_dens !=0 ),
                    x ="mean_icseq",
                    y = "mean_icfrust",
                    color = "sd_dens",
                    size = "n_subgenus",
                    shape = 'cl_phylo',
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"),
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "pearson", label.x = 2.5,label.y = 1, label.sep = "\n"),
                    xlab = "SeqIC",
                    ylab= "FrustIC"
)+
  theme_minimal()+
  theme(legend.position = "right")+
  scale_color_viridis(discrete = FALSE, option = "D",breaks = c(0.1, 0.6),
                      labels = c("0.1 (Balanced)","0.6 (Imbalanced)"))+
  scale_shape_manual(values = c(15,17,16,8))+
  guides(color = guide_colourbar(title = "Phylogenetic\ndiversity",ncol=2),
         size = guide_legend(title = "# subgenus species",ncol=2), nrow = 2,
         shape = guide_legend(title = "Genus\nclassification",ncol=1))#+
  


scatter

saveRDS(scatter, 'figureS8C.RDS')



