library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
library(bio3d)
library(tidyr)
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



nseqs_plddt = fread('../../figureS6/nseqs_plddt.txt')
select = subset(nseqs_plddt, nseqs >= 10 & y >= 80)
colnames(select)[colnames(select) == "cluster"] = "Cluster"

res_df = fread('../../figureS7/frust_3dmapper_sdps_results.txt')
colnames(res_df)[colnames(res_df) == "Gene"] = "prot"


res_df_filtered = merge(res_df, select)

res_df_filtered_sum = subset(res_df_filtered, subcluster == "rest") %>% group_by(prot, Cluster, PDB_chain_length) %>%
  summarise(mean_icseq = mean(ICseq),
            mean_icfrust = mean(ICTot))


list_plddt = list()
list_iupred = list()
for (p in unique(res_df_filtered_sum$prot)) {
  print(p)
  
  plddt_file = Sys.glob(paste(wd,
                              '/',
                              p,
                              "/",
                              p,
                              "_mean_plddt_plpro_models.txt",
                              sep = ""))
  
  plddt = fread(plddt_file)
  plddt$prot = p
  colnames(plddt)[colnames(plddt) == "ali_pos"] = "Subalignment_position"
  colnames(plddt)[colnames(plddt) == "cluster"] = "Cluster"
  
  plddt$n = 1
  plddt_sum = plddt %>% group_by(prot, Cluster) %>% filter(mean_plddt < 80) %>% summarise(n_lowplddt = n())
  iupred_sum = plddt %>% group_by(prot, Cluster) %>% filter(mean_iupred > 0.5) %>% summarise(n_highiupred = n())
  
  list_plddt[[p]] = plddt_sum
  list_iupred[[p]] = iupred_sum
  
}


model_quality = rbindlist(list_plddt)
flexible_regions = rbindlist(list_iupred)

res_df_sum_filtered_model_quality = left_join(res_df_filtered_sum, model_quality, by =
                                                c("prot", "Cluster"))
res_df_sum_filtered_model_quality$dens_lowplddt = res_df_sum_filtered_model_quality$n_lowplddt / res_df_sum_filtered_model_quality$PDB_chain_length

res_df_sum_filtered_flexible_regions = left_join(res_df_filtered_sum, flexible_regions, by =
                                                   c("prot", "Cluster"))
res_df_sum_filtered_flexible_regions$dens_highiupred = res_df_sum_filtered_flexible_regions$n_highiupred / res_df_sum_filtered_flexible_regions$PDB_chain_length


scatter = ggscatter(
  res_df_sum_filtered_model_quality,
  x = "mean_icseq",
  y = "mean_icfrust",
  color = "dens_lowplddt",
  add = "reg.line",
  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE,
  # Add confidence interval
  cor.coef = TRUE,
  # Add correlation coefficient. see ?stat_cor
  cor.coeff.args = list(
    method = "pearson",
    label.x = 2.5,
    label.y = 1,
    label.sep = "\n"
  ),
  xlab = "SeqIC",
  ylab = "FrustIC"
) +
  theme_minimal() +
  theme(legend.position = "top")  +
  scale_color_viridis(discrete = FALSE,
                      option = "mako",
                      direction = -1) +
  guides(color = guide_colourbar(title = "Proportion residues\n with pLDDT < 80"))#+


scatter
saveRDS(scatter, 'figure8A.RDS')

