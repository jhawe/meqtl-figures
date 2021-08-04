# ------------------------------------------------------------------------------
# Plots figure 1A - CHessboard overview plot of cosmo pairs
# Needs the prepare_chessboard_data.R script to be run beforehand
# ------------------------------------------------------------------------------
require(plotrix)

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")

# load prepared plot data
load('chessboard/preprocessed_cis_rounded_100k.RData')

# plot -------------------------------------------------------------------------

pdf(file = "chessboard.pdf",
    height = 18,
    width = 18)

par(
  mar = c(5,4,4,4),
  mgp = c(2.5, 1, 0),
  oma = c(5, 5, 0, 0),
  cex.axis = 1.5,
  col.axis = 'dark grey'
)

mat <- matrix(nrow = 2, ncol = 2)
mat[1:2, 1:2] <- seq(1, 4, 1)
layout(mat, widths = c(0.5, 3), heights = c(3, 0.5))
lwid <- 2.5

#SNP density
plot(
  c(0, max(res.dens$sn)),
  c(chrborders[1], chrborders[23]),
  col = 'white',
  xlab = '',
  ylab = '',
  yaxt = 'n',
  xaxt = 'n',
  yaxs = 'i'
)

for (chr in 1:22) {
  dens = res.dens[res.dens$c == chr, c('g', 'sn')]
  points(
    dens$sn,
    dens$g,
    type = 'l',
    col = c('black', 'grey')[(chr %% 2) + 1],
    lwd = lwid
  )
  abline(
    h = chrborders[chr],
    col = 'dark grey',
    lty = 2 ,
    lwd = lwid
  )
}

axis(2, tckpos, seq(1, 22, 1), tick = FALSE, las = 1)
axis(2, chrborders, rep('', 23))
axis(3, c(0, 20000), c(0, '20k'))
mtext(
  'SNPs',
  side = 2,
  line = 0,
  adj = 0.55,
  cex = 3,
  outer = T
)
plot(
  0,
  xaxt = 'n',
  yaxt = 'n',
  bty = 'n',
  pch = '',
  ylab = '',
  xlab = ''
)


#main plot
plot(
  c(chrborders[1], chrborders[23]),
  c(chrborders[1], chrborders[23]),
  xlab = '',
  ylab = '',
  xaxt = 'n',
  yaxt = 'n',
  col = 'white',
  xaxs = 'i',
  yaxs = 'i'
)

for (m in 1:(length(chrborders) - 1)) {
  i = ifelse(m %% 2 == 0, 1, 2)
  for (s in 1:(length(chrborders) - 1)) {
    rect(
      chrborders[m],
      chrborders[s],
      chrborders[m + 1],
      chrborders[s + 1],
      col = color[i],
      border = NA
    )
    i = i + 1
  }
}
print(date())
points(sigtab.cis[, 1],
       sigtab.cis[, 2],
       pch = 20,
       col = "#009E73",
       cex = 1.2)
points(
  sigtab.longrange[, 1],
  sigtab.longrange[, 2],
  pch = 20,
  col = "#AA4499",
  cex = 1.2
)
points(
  sigtab.trans[, 1],
  sigtab.trans[, 2],
  pch = 20,
  col = 'black',
  cex = 1.2
)
print(date())
box(lwd = 1.5)

#CpG density
plot(
  c(chrborders[1], chrborders[23]),
  c(0, max(res.dens$mn)),
  col = 'white',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  yaxt = 'n',
  xaxs = 'i'
)
for (chr in 1:22) {
  dens = res.dens[res.dens$c == chr, c('g', 'mn')]
  points(
    dens$g,
    dens$mn,
    type = 'l',
    col = c('black', 'grey')[(chr %% 2) + 1],
    lwd = lwid
  )
  abline(
    v = chrborders[chr],
    col = 'dark grey',
    lty = 2 ,
    lwd = lwid
  )
}
axis(1, tckpos, seq(1, 22, 1), tick = FALSE, padj = 1)
axis(1, chrborders, rep('', 23))
axis(2, c(0, 3000), c(0, '3k'), las = 1)
mtext(
  'CpGs',
  side = 1,
  line = 0,
  adj = 0.5,
  cex = 3,
  outer =T
)

dev.off()
