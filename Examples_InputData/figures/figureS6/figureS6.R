library(data.table)
library(ggplot2)
library(ggnewscale)
library(seqinr)
library(peRReo)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(yarrr)
library(RColorBrewer)
library(data.table)
library(stringr)
library(dplyr)
library(plyr)
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

########
# pLDDT plot
########
plddt_fs = Sys.glob(file.path(wd, "*/*_mean_plddt_plpro_models.txt"))

list_results = list()
for (f in plddt_fs) {
  protname = basename(dirname(f))
  plddt = fread(f)
  plddt$prot = protname
  list_results [[protname]] = plddt
}

df = rbindlist(list_results)
df_sd = df %>% group_by(cluster, prot) %>% summarise(y = mean(mean_plddt, na.rm = T),
                                                     sd = sd(mean_plddt, na.rm = T))

df_sd = df_sd %>% mutate(
  ymin = y - sd,
  ymax = y + sd,
  sdp_cluster = as.character(cluster)
)



df_sd = left_join(df_sd, prot_df)

df_sd$protein = factor(df_sd$protein, levels =  rev(
  c(
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
    'E',
    'N_Cterm',
    'N_Nterm',
    'S_RBD',
    'orf3a',
    'orf7a',
    'orf8',
    'orf9b'
  )
))


p = ggplot(
  subset(df_sd, protein != "S_RBD"),
  aes(
    x = protein,
    y = y,
    ymin = ymin,
    ymax = ymax,
    color = sdp_cluster
  )
) +
  geom_pointrange() +
  theme_minimal() +
  xlab("") +
  ylab("mean pLDDT") +
  coord_flip() +
  geom_hline(yintercept = 80, linetype = "dashed") +
  scale_color_manual(values = brewer.pal(7, 'Paired'))


#######################
# NUMBER OF SEQS PLOT
#######################


dir = Sys.glob(file.path(wd, "*"))
results = list()
results2 = list()
for (i in 1:length(dir)) {
  trimmed = Sys.glob(file.path(dir[i], 'trimmed*'))
  clusters = Sys.glob(file.path(dir[i], '*_cluster*txt'))
  if (length(trimmed) > 0 & length(clusters) > 0) {
    cl = fread(clusters)
    cl$seq_name = str_replace_all(cl$seq_name, '\\/', '_')
    colnames(cl)[1] = "protein"
    dt = fread(trimmed)
    prot =  basename(dirname(trimmed))
    results[[i]] = data.frame(nseqs = nrow(dt), prot)
    dt2 = dt
    dt2$prot = prot
    dt2 = left_join(dt2, cl)
    results2[[i]] = dt2
  }
}

r2 = rbindlist(results2)
r2$n = 1
r3 = r2 %>% dplyr::group_by(prot, cluster) %>%
  dplyr::summarise(nseqs = sum(n))

ts = merge(df_sd, r3)


r3$cluster = factor(r3$cluster, levels = c("1", "2", "3", "4", "5", "6", "7"))
r3 = r3[order(r3$cluster), ]
df_cumsum = r3 %>% group_by(prot) %>% arrange(cluster)   %>% dplyr::mutate(label_ypos = cumsum(nseqs))

df_cumsum = left_join(df_cumsum, prot_df)

df_cumsum$protein = factor(df_cumsum$protein, levels =  rev(
  c(
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
    'S_RBD',
    "S_protein",
    'E',
    'N_Cterm',
    'N_Nterm',
    'orf3a',
    'orf7a',
    'orf8',
    'orf9b'
  )
))



df_cumsum$label2 = ifelse(df_cumsum$nseqs < 11, df_cumsum$nseqs, NA)
df_cumsum$cluster = factor(df_cumsum$cluster, levels = rev(c("1", "2", "3", "4", "5", "6", "7")))

g1 = ggplot(subset(df_cumsum, !is.na(cluster)),
            aes(
              y = protein,
              x = nseqs,
              fill = cluster,
              label = nseqs
            )) +
  geom_bar(stat = "identity") +                                                              # Add values on top of bars
  geom_text(size = 3.5,
            position = position_stack(vjust = 0.5),
            color = "white") +
  theme_minimal() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(
    title = "# SDP cluster",
    nrow = 1,
    reverse = TRUE
  )) +
  ylab("SARS-CoV-2 protein") +
  xlab("Number of sequences") +
  scale_fill_manual(values =  rev(brewer.pal(7, 'Paired')))#+
legend = get_legend(g1)



c = plot_grid(
  legend,
  plot_grid(
    g1 + theme(legend.position = "None"),
    p + theme(legend.position = "None"),
    ncol = 2,
    rel_widths = c(2, 1),
    labels = c("A", "B")
  ),
  ncol = 1,
  rel_heights = c(0.1, 1.5)
)

svg('figureS6.svg', width = 12, height = 8)
c
dev.off()
