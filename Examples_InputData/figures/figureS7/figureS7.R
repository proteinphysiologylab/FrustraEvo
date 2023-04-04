library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ggpubr)
##################################

wd = "../figure5/SARS2_FrustraEvo_input_and_output"

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

select = c(
  'NSP2_protein',
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
  "S_protein"
)



nseqs_plddt = fread('../figure5/nseqs_plddt.txt')
select = subset(nseqs_plddt, nseqs >= 10 & y >= 80)
colnames(select)[colnames(select) == "cluster"] = "Cluster"


########
# SDPs vs Frust
########
mapper_fs = Sys.glob(file.path(wd, "*/*mapped_positions_unique.csv"))

# read 3dmapper results and merge them with frustration single res results
frust_3dmapper_results = list()
for (f in mapper_fs) {
  prot = basename(dirname(dirname(dirname(f))))
  frustraevo_fs = Sys.glob(
    file.path(
      wd,
      prot,
      "frustraevo_output/*/OutPut*/Equivalences/IC_seqfrust.csv"
    )
  )
  
  m = fread(f)
  m_ref = subset(m, IsRef == "yes" & Pident == 100)
  
  cluster_res = list()
  for (f2 in frustraevo_fs) {
    cluster = str_remove(basename(dirname(dirname(dirname(
      f2
    )))), 'cluster')
    if (cluster == "all") {
      subcluster = "all"
      if (prot != "S_protein" & prot != "S_RBD") {
        cluster = 1
      } else {
        cluster = 4
        
      }
    } else {
      cluster = as.numeric(cluster)
      subcluster = "rest"
    }
    ff = fread(f2)
    ff$Cluster = cluster
    ff$subcluster = subcluster
    colnames(ff)[colnames(ff) == "Res"] = "Trimmed_protein_position"
    cluster_res[[paste(cluster, subcluster, sep = "_")]] = ff
  }
  
  frustra_all = rbindlist(cluster_res)
  
  df = unique(merge(m_ref, frustra_all))
  frust_3dmapper_results[[prot]] = df
}


res_df = rbindlist(frust_3dmapper_results)
fwrite(res_df, 'frust_3dmapper_sdps_results.txt')
#######################################
# FINAL FIGURE
#######################################
colnames(res_df)[colnames(res_df) == "Gene"] = "prot"
res_df_sum = subset(res_df, subcluster == "rest") %>% group_by(prot, Cluster) %>%
  summarise(mean_icseq = mean(ICseq),
            mean_icfrust = mean(ICTot))

res_df_sum_filtered = merge(res_df_sum, select)
res_df_sum_filtered$lab = paste(res_df_sum_filtered$protein,
                                res_df_sum_filtered$Cluster,
                                sep = "_")


scatter = ggscatter(
  subset(res_df_sum_filtered, prot != "S_RBD"),
  x = "mean_icseq",
  y = "mean_icfrust",
  color = "protein",
  add = "reg.line",
  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE,
  # Add confidence interval
  cor.coef = TRUE,
  # Add correlation coefficient. see ?stat_cor
  cor.coeff.args = list(
    method = "pearson",
    label.x = 2.25,
    label.y = 1.15,
    label.sep = "\n"
  ),
  ylim = c(0.2, 1.2),
  xlab = "SeqIC",
  ylab = "FrustIC"
) +
  theme_minimal() +
  geom_label_repel(
    aes(label = lab, color = protein),
    size = 4,
    max.overlaps = 1000,
    box.padding = 2
  ) +
  theme(
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(0.1, "cm"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


svg('figureS7.svg', width = 11, height = 10)
scatter
dev.off()
