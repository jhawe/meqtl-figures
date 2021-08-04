# create main figure 6 from individual figures

library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
theme_set(theme_cowplot())

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")

# cannot be loaded in R4+ -> save in R3X then manually copy in
#load("iqtl_figure/asthma_example.RData")
#asthma.example

# plot saved in 'cosmo.repl'
load("data/iqtl_figure/cosmo_replication_isolated_cells.RData")
cosmo.repl <- cosmo.repl + theme(legend.position = c(0.75,0.8))
cosmo.repl

# plot saved in "manhatten"
load("data/iqtl_figure/example_iqtl_MAGNETIC_HDL.C_manhattan.RData")
manhattan <- manhattan + background_grid() + 
  theme(panel.border = element_rect(size = 1, color="black"))

# plot saved in "enrich.plot"
load("data/iqtl_figure/volcano_plot.RData")
enrich.plot

# iQTL replication
load("data/iqtl_figure/cohort_replication_cd8t.RData")
rm(all.iqtl)
cohort_replication <- cohort_replication + 
  theme(legend.position = "none")

# US legal
pdf("figure6.pdf", width=14, height=8.5)
ggarrange(ggarrange(cohort_replication, cosmo.repl, enrich.plot, 
                    ncol=3, widths = c(1.5,2.5,2.5), labels = "AUTO"),
          ggarrange(NULL, manhattan, ncol=2, labels=c("D", "E")),
          nrow=2)
dev.off()
