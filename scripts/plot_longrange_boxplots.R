#
# Plot longrange boxplots
#

library(tidyverse)
library(cowplot)
library(scales)

theme_set(theme_cowplot())

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")

# ------------------------------------------------------------------------------
# boxplots for cpg enrichment in distance groups
load(paste('data/LONGRANGE/', 1, '_v2.RData', sep=''))
ratios=as.numeric(as.character(resu[,3]))
for(i in 2:1100){
  tryCatch({
    load(paste('data/LONGRANGE/', i, '_v2.RData', sep=''))
    ratios=cbind(ratios, as.numeric(as.character(resu[,3])))
  }, error = function(error) {return(NA)})
}
ratios=ratios[,1:1000]

#plot
resu=rep(NA, (nrow(ratios)-1))
ratios=ratios[c(1:10,13,12,11),]
for(i in 1:(nrow(ratios)-1)){
  resu[i]=t.test(ratios[i,],ratios[13,], alternative ='greater')$p.value
}

pdf("longrange_distances_cpg_enrichment.pdf", width = 17, height=8)

par(mar = c(4, 8, 2, 2), mgp = c(4.5, 1, 0))
boxplot(
  t(ratios + 1E-10),
  log = 'y',
  col = c(rep('lightblue', 10), rep('orange', 2), 'blue'),
  xaxt = 'n',
  xlab = '',
  ylab = 'Proportion of Significant CpGs',
  cex.lab = 2,
  cex.axis = 1.6,
  ylim = c(1E-10, 1E-1)
)
axis(
  1,
  1:13,
  c(
    '<1Mb',
    '1-2Mb',
    '2-3Mb',
    '3-4Mb',
    '4-5Mb',
    '5-6Mb',
    '6-7Mb',
    '7-8Mb',
    '8-9Mb',
    '9-10Mb',
    '>1Mb',
    '>10Mb',
    'trans'
  ),
  cex.axis = 1.5
)
b = boxplot(t(ratios + 1E-10), log = 'y', plot = F)
text(
  1:12,
  b$stats[5, 1:12],
  ifelse(
    resu > 0,
    format(
      resu,
      scientific = T,
      nsmall = 2,
      digits = 1
    ),
    '<1e-324'
  ),
  col = 'black',
  pos = 3,
  offset = 2,
  cex = 1.3,
  font = 2
)

dev.off()
