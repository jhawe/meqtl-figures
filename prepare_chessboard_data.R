# Preprocess data to be used in chessboard plotting (also reduce RData sizes)

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")

# load all raw/previous data from BL and ML
load('chessboard/sigproces_110214.RData')
load('chessboard/anno_170714.RData')
load("chessboard/cosmopairs_chessboard_070617.RData")

rm(sig.mat, stats, win.mat, dens)
rm(annometh, annosnp, chrend, chrstart, cumuchrl2, xcld)

gc(full=T)

# prepare plotting data --------------------------------------------------------
chrborders = cumuchrl[1:23]
tckpos = rep(NA, 22)
tckpos[1] = chrborders[2] / 2
for (i in 2:22) {
  tckpos[i] = chrborders[i] + (chrl[i] / 2)
}

# colors for background grid
color <- rep(c('#0000FF20', 'white'), (length(chrborders) + 1))
mn <- win$mn
sn <- win$sn
win.dens <- cbind(mn, sn)
res.dens <- data.frame(cbind(win$c,
                             ceiling(((win$egw - win$sgw) / 2
                             ) + win$sgw),
                             win.dens))
colnames(res.dens)[1:2] = c('c', 'g')

save.image("chessboard/preprocessed.RData")

# cis and longrange cis associations are too many to properly display in a
# meaningful way -> round positions and save unique sets

sigtab.cis$gpos.cpg <- round(sigtab.cis$gpos.cpg,-3)
sigtab.cis$gpos.snp <- round(sigtab.cis$gpos.snp,-3)
sigtab.cis <- unique(sigtab.cis)

sigtab.longrange$gpos.cpg <- round(sigtab.longrange$gpos.cpg,-3)
sigtab.longrange$gpos.snp <- round(sigtab.longrange$gpos.snp,-3)
sigtab.longrange <- unique(sigtab.longrange)

save.image("chessboard/preprocessed_cis_rounded_1k.RData")

sigtab.cis$gpos.cpg <- round(sigtab.cis$gpos.cpg,-4)
sigtab.cis$gpos.snp <- round(sigtab.cis$gpos.snp,-4)
sigtab.cis <- unique(sigtab.cis)

sigtab.longrange$gpos.cpg <- round(sigtab.longrange$gpos.cpg,-4)
sigtab.longrange$gpos.snp <- round(sigtab.longrange$gpos.snp,-4)
sigtab.longrange <- unique(sigtab.longrange)

save.image("chessboard/preprocessed_cis_rounded_10k.RData")


sigtab.cis$gpos.cpg <- round(sigtab.cis$gpos.cpg,-5)
sigtab.cis$gpos.snp <- round(sigtab.cis$gpos.snp,-5)
sigtab.cis <- unique(sigtab.cis)

sigtab.longrange$gpos.cpg <- round(sigtab.longrange$gpos.cpg,-5)
sigtab.longrange$gpos.snp <- round(sigtab.longrange$gpos.snp,-5)
sigtab.longrange <- unique(sigtab.longrange)

save.image("chessboard/preprocessed_cis_rounded_100k.RData")

# all done
sessionInfo()