## LOAD PACKAGES
warn.conflicts = FALSE
options(echo = FALSE, verbose = F,warn = -1) 

#also required Biostrings!!!
requiredPackages = c('ggplot2', 'ggnewscale', 'cowplot',
                     'seqinr','stringr', 'data.table', "dplyr", "tidyr","ggdendro","ggtext",'GetoptLong') 
suppressMessages(
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE))
      install.packages(p)
    library(p, character.only = TRUE)
  }
)
