library(cowplot)
library(ggplot2)

list_corr_pdz = lapply(Sys.glob("../figure2/PDZ3/PSD95_corr_ic*_vs_*fitness.RDS"), readRDS)

plot_sh3_logo = lapply(c("../figure2/SH3/bplot.RDS",
                         "../figure2/SH3/legend_bplot.RDS",
                         "../figure2/SH3/legend_logo.RDS",
                         "../figure2/SH3/slogo.RDS"), readRDS)

list_corr_sh3 = lapply(Sys.glob("../figure2/SH3/SH3_corr_ic*_vs_*fitness.RDS"), readRDS)

list_corr_kras = lapply(Sys.glob("../figure2/KRAS/RAS_Rojas_et_al/KRAS_corr_ic*_vs_*fitness.RDS"), readRDS)

list_corr_icfrust_icseq = lapply(Sys.glob(c("../figure2/*/*_corr_ic*_vs_*ic*.RDS",
                                          "../figure2/*/*Rojas*/*_corr_ic*_vs_*ic*.RDS")), readRDS)

library(cowplot)
svg('figureS2.svg', width = 11, height = 11)
plot_grid(list_corr_sh3[[2]]+ labs(color ="Frust. State"),
          list_corr_sh3[[4]]+ labs(color ="Frust. State"),
          list_corr_icfrust_icseq[[2]]+ labs(color ="Frust. State"),
          list_corr_pdz[[2]]+ labs(color ="Frust. State"),
          list_corr_pdz[[4]]+ labs(color ="Frust. State"),
          list_corr_icfrust_icseq[[1]]+ labs(color ="Frust. State"),
          list_corr_kras[[2]]+ labs(color ="Frust. State"),
          list_corr_kras[[4]]+ labs(color ="Frust. State"),
          list_corr_icfrust_icseq[[3]]+ labs(color ="Frust. State"),
          align = 'hv',
          labels = c("A","B","C", "D", "E","F","G", "H", "I"),
          label_size = 20,
          ncol = 3
)
dev.off()
