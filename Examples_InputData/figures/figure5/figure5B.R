library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
library(ggh4x)
##################################
wd = "../SARS2_FrustraEvo_input_and_output"

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
    "S_protein",
    'S_RBD',
    'E',
    'N_Cterm',
    'N_Nterm',
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
    "S_protein",
    'S_RBD',
    'E_protein',
    'N_C-terminaldomain',
    'N_N-terminaldomain',
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

nseqs_plddt = fread('../nseqs_plddt.txt')
select = subset(nseqs_plddt, nseqs >= 10 & y >= 80)
colnames(select)[colnames(select) == "cluster"] = "Cluster"
########
# SDPs vs Frust
########
mapper_fs =  Sys.glob(file.path(wd, "*/*mapped_positions_unique.csv"))

# read 3dmapper results and merge them with frustration single res results
frust_3dmapper_results = list()
for (f in mapper_fs) {
  prot = basename(dirname(f))
  frustraevo_fs = Sys.glob(
    file.path(
      wd,
      prot,
      "frustraevo_output/*/OutPut*/Equivalences/IC_seqfrust.csv"
    )
  )
  
  m = fread(f)
  m_ref = subset(m, IsRef == "yes" & Pident == 100)
  if (prot == "S_protein" | prot == "S_RBD") {
    m_ref = m_ref %>% group_by(Protein_position, Cluster) %>% filter(Trimmed_protein_position == min(Trimmed_protein_position))
  }
  cluster_res = list()
  for (f2 in frustraevo_fs) {
    cluster = as.numeric(str_split_fixed(basename(dirname(dirname(
      f2
    ))), 'cluster', 2)[, 2])
    ff = fread(f2)
    ff$Cluster = cluster
    colnames(ff)[colnames(ff) == "Res"] = "Trimmed_protein_position"
    cluster_res[[cluster]] = ff
  }
  
  frustra_all = rbindlist(cluster_res)
  
  df = unique(merge(m_ref, frustra_all))
  frust_3dmapper_results[[prot]] = df
}


frust_3dmapper = rbindlist(frust_3dmapper_results)


# subset for cluster containing sars2 sequence
sarbecov = subset(frust_3dmapper,
                  ((Cluster == 1 & Gene != "S_protein") |
                     (
                       Cluster == 4 &
                         Structure_feature_id == "SARS-CoV2.pdb__oldS_SARS-CoV2_A_NA_NA"
                     )
                  ) &
                    Gene %in% subset(select, Cluster == 1)$prot)

# fix s protein pdb length
sarbecov$PDB_chain_length = ifelse(sarbecov$Gene == "S_protein", 1131 , sarbecov$PDB_chain_length)

# add labels to sdps and ICfrust
#sarbecov$sdp = ifelse(sarbecov$Consequence =="1st_level_SDP", "yes", "no" )
sarbecov$ICFrust = ifelse(sarbecov$ICTot > 0.5, sarbecov$FrustEstado, 'ICFrust≤0.5')

# select only proteins that passed the plddt and nseqs filter
colnames(sarbecov)[colnames(sarbecov) == "Gene"] = "prot"
sarbecov = merge(sarbecov, prot_df)

sarbecov_sum = sarbecov %>% group_by(protein, ICFrust, PDB_chain_length) %>%
  summarise(count = n())

sarbecov_sum$density = sarbecov_sum$count / sarbecov_sum$PDB_chain_length


levs = (
  subset(sarbecov_sum, ICFrust == 'ICFrust≤0.5') %>%
    group_by(protein) %>%
    summarise(sum_dens = sum(density)) %>%
    arrange(desc(sum_dens))
)$protein

sarbecov_sum$protein = factor(sarbecov_sum$protein, levels = levs)
sarbecov_sum$ICFrust = factor(sarbecov_sum$ICFrust, levels = rev(c('ICFrust≤0.5', "NEU", "MIN", "MAX")))



p1 = ggplot(sarbecov_sum, aes(y = protein, x = density, fill = ICFrust)) +
  geom_bar(stat = "identity",
           color = "black")   +
  scale_fill_manual(
    values = c(
      'ICFrust≤0.5' = "lightblue",
      "MAX" = "red",
      "MIN" = "green3",
      "NEU" = "grey70"
    ),
    breaks = rev(c('ICFrust≤0.5', "MAX", "MIN", "NEU"))
  ) +
  scale_alpha_manual(values = c(.5, 1)) +
  theme_grey() +
  xlab("Density") +
  ylab("SARS-CoV-2 Protein") +
  labs(fill = "ICFrust > 0.5") +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_x_continuous(expand = c(0, 0)) +
  force_panelsizes(rows = c(1, 2))

saveRDS(p1, 'figure5B.RDS')
