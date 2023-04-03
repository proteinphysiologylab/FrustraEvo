library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(seqinr)
library(stringr)
library(dplyr)
##################################
corrplot = readRDS("./figure5A/figure5A.RDS")
barplot = readRDS("./figure5B/figure5B.RDS")
msfa = readRDS("./figure5C/figure5C.RDS")


corrplot = corrplot +
  theme_classic()+
  theme(legend.position = c(0.54, 0.17),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.key.size = unit(0.1, "cm"),
        text = element_text(size = 12) )+
  ylab("FrustIC")+
  xlab("SeqIC")

barplot = barplot+
  labs(fill="FrustIC>0.5")+
  scale_fill_manual(values=c("red","green3", "grey70",  "lightblue"),
                      labels = c("MAX","MIN","NEU", "FrustIC<=0.5"))


msfa = msfa+
  labs(size ="SeqIC")+
  scale_fill_gradientn(values = c(1, .79, .5, 0.11,0.0909, 0), 
                       colours=c("green", "grey", "grey", "grey","red", "red"),
                       name= "SRFI median index",
                       limits=c(-1.2, 1.1),
                       oob = scales::squish,
                       na.value = "white")

svg("figure5.svg", width = 13, height = 8 )
plot_grid(plot_grid(
                    corrplot,
                    align = "hv",
          barplot,
          ncol = 2, 
          rel_widths = c(.75,1),
          labels = c("A", "B"), label_size = 20),
          NULL,
          msfa,
          align = "hv",
          nrow =3,
          rel_heights = c(0.95,0.05,0.5),
          labels = c("","", "C    PLPro"), label_size = 20)

dev.off()
