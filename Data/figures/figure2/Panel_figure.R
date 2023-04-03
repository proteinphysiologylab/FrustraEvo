library(cowplot)
library(ggplot2)

list_corr_pdz = lapply(Sys.glob("./PDZ3/PSD95_corr_ic*_vs_*fitness.RDS"), readRDS)

plot_sh3_logo = lapply(c("./SH3/bplot.RDS",
                         "./SH3/legend_bplot.RDS",
                         "./SH3/legend_logo.RDS",
                         "./SH3/slogo.RDS"), readRDS)

list_corr_sh3 = lapply(Sys.glob("./SH3/SH3_corr_ic*_vs_*fitness.RDS"), readRDS)

list_corr_kras = lapply(Sys.glob("./KRAS/RAS_Rojas_et_al/KRAS_corr_ic*_vs_*fitness.RDS"), readRDS)


######
# FINAL FIGURES
#######

svg('figure2.svg', width = 17, height = 8)
plot_grid(
  #NULL,
  plot_grid(
    NULL ,
    plot_grid(NULL, labels = c("A"),label_size = 20),
   # plot_grid(NULL,
              plot_grid(
                plot_grid(NULL,plot_sh3_logo[[3]], NULL, ncol = 3) ,
                plot_sh3_logo[[4]]+ylab("SeqIC"),
                plot_grid(NULL,plot_sh3_logo[[2]],NULL, ncol = 3) ,
                plot_sh3_logo[[1]]+ylab("FrustIC")+theme(legend.position = "none"),
                align= "v",axis = "lr",
                rel_heights = c(0.1,0.75,0.05,1),
                ncol = 1),
    #          rel_widths = c(0.08,1),
    #          ncol = 2),
    plot_grid(list_corr_sh3[[3]]+ labs(color ="Frustration State"),
              list_corr_sh3[[1]]+ labs(color ="Frustration State"),
              NULL,
              rel_widths = c(1,1,0.05),
              nrow=1,
              align = 'hv', labels = c("B","C"),label_size = 20),
    nrow = 4,
    labels =c("                       GRB2-SH3"),
    label_size = c(23),
    rel_heights = c(0.1,0.05,1,1))+
    theme(panel.border = element_rect(color = "black", fill=NA, size  =1)),
  #NULL,
  plot_grid(
    plot_grid(NULL,
              plot_grid(
                list_corr_pdz[[3]]+ labs(color ="Frustration State"),
                list_corr_pdz[[1]]+ labs(color ="Frustration State"),
                NULL,
                ncol = 3,
                rel_widths = c(1,1,0.05),
                
                align = 'hv',
                labels = c("D", "E"),label_size = 20),
              ncol = 1,
              rel_heights = c(0.1,1),
              labels =c("                  PSD95-PDZ3"),
              label_size = c(23))+
      theme(panel.border = element_rect(color = "black", fill=NA, size  =1)),
    
    plot_grid(NULL, 
              plot_grid(list_corr_kras[[3]]+ labs(color ="Frustration State"),
                        list_corr_kras[[1]]+ labs(color ="Frustration State"),
                        NULL,
                        ncol = 3,
                        rel_widths = c(1,1,0.05),
                        align = 'hv',
                        labels = c("F", "G"),label_size = 20
              ),
              ncol = 1,
              rel_heights = c(0.1,1),
              labels =c("                          KRAS"),
              label_size = c(23))+
      theme(panel.border = element_rect(color = "black", fill=NA, size  =1)),
    ncol = 1
  ), 
  label_size = 14,
  ncol = 2,#4,
  rel_widths = c(1,1)
)
dev.off()


