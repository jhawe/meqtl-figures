#
# PLot the cis histogram for EDF2
#

library(tidyverse)
library(cowplot)
library(scales)

theme_set(theme_cowplot())

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")

load("data/cosmopairs_combined_151216.RData")

# ------------------------------------------------------------------------------
# histogram of significant associations
require(plotrix)

sigtab.cis=cosmo[cosmo$snp.chr==cosmo$cpg.chr & abs(cosmo$snp.pos-cosmo$cpg.pos)<1000000,]
sigtab.cis$distance=sigtab.cis$snp.pos-sigtab.cis$cpg.pos
tmp=sigtab.cis[with(sigtab.cis, order(cpg, abs(distance))),]
sigtab.cis=tmp[!duplicated(tmp$cpg),]

pdf("cishist.pdf", width = 13)
par(
  mar = c(7, 7, 2, 2),
  mgp = c(4.2, 0.8, 0),
  cex.lab = 2,
  cex.axis = 1.7
)
h = hist(sigtab.cis$distance,
         breaks = seq(-1010000, 1010000, 20000),
         plot = F)
colors = rep('lightgrey', length(h$counts) + 2)
gap.barplot(
  c(0, h$counts + 1, 0),
  gap = c(1800, 60500),
  xlab = 'Distance SNP-CpG (kb)',
  col = colors,
  main = '',
  ylab = '# CpGs with Significant Associations',
  xaxt = 'n',
  yaxt = 'n',
  ylim = c(0, 2800),
  cex.lab = 2
)
abline(h = 1870, col = 'white', lwd = 35)
box(fg = 'grey')
axis.break(2, 1900, brw = 0.05, style = 'zigzag')
axis(1, seq(2, 102, 10), gsub(' ', '', format(seq(-1000, 1000, 200), big.mark =
                                                ',')), fg = 'grey')
axis(2, seq(0, 1500, 500),  cex.axis = 1.7)
axis(2, c(2300, 2800), c(61000, 61500), cex.axis = 1.7)
dev.off()
