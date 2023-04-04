#plot_pdz_logo  = readRDS('/media/victoria/VICTORIA_EXTERNAL_DISK/CORONAVIRUS/frustration/paper_figs/lehner/RDS/plot_logo_PDZ.RDS')
library(ggplot2)
plddt= readRDS("./figureS8A/figureS8A.RDS")
iupred= readRDS("./figureS8B/figureS8B.RDS")
density = readRDS("./figureS8C/figureS8C.RDS")

library(cowplot)
######
# FINAL FIGURE
#######

svg('figureS8.svg', width = 10, height = 8)
plot_grid(
  plot_grid(plddt,
            iupred,
            ncol =2,
            labels = c("A","B"),
            label_size = 20,
            align = "h"),
  density,
  labels = c("","C"),
  label_size = 20,
  ncol = 1,#4,
  rel_heights =  c(0.75,1),
  align = "v"
)
dev.off()


