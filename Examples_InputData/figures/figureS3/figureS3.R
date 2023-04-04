
list_corr_kras = lapply(Sys.glob("../KRAS/KRAS_corr_ic*_vs_*fitness.RDS"), readRDS)
list_corr_icfrust_icseq = lapply(Sys.glob("../KRAS/*_corr_ic*_vs_*ic*.RDS"), readRDS)

library(cowplot)
library(ggplot2)

#####
# FINAL FIGURES
#######

svg('figureS3_missing.svg', width = 12, height = 8)
plot_grid(
  list_corr_kras[[1]],
  list_corr_kras[[3]],
  NULL,
  list_corr_kras[[2]],
  list_corr_kras[[4]],
  list_corr_icfrust_icseq[[1]],
  labels= c("B", "C","D",
          "E", "F", "G"),
  ncol=3,align = "vh",
  label_size = 22
)
  
dev.off()
  
  