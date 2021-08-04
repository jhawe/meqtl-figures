require(plotrix)

#get data

#contains marker densities
load('chessboard/sigproces_110214.RData')  
load('TODO')

load('cosmopairs_combined_151216.RData')
chrborders = cumuchrl[1:23]
tckpos = rep(NA, 22)
tckpos[1] = chrborders[2] / 2
for (i in 2:22) {
  tckpos[i] = chrborders[i] + (chrl[i] / 2)
}
color <- rep(c('#0000FF20', 'white'), (length(chrborders) + 1))
mn <- win$mn
sn <- win$sn
win.dens <- cbind(mn, sn)
res.dens <- data.frame(cbind(win$c, 
                             ceiling(((win$egw - win$sgw) / 2) + win$sgw),
                             win.dens))
colnames(res.dens)[1:2] = c('c', 'g')


# plot -------------------------------------------------------------------------

pdf(file="/home/mcloh/chessboard.pdf", height=18, width=18)

par(mar=c(3,3,3,3), 
    mgp = c(2.5,1,0), 
    oma=c(30,30,10,10), 
    cex.axis=5, 
    col.axis='dark grey', 
    cex.lab=6)

mat <- matrix(nrow = 2, ncol = 2)
mat[1:2, 1:2] <- seq(1, 4, 1)
layout(mat, widths = c(0.3, 3), heights = c(3, 0.3))
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

for(chr in 1:22){
  dens=res.dens[res.dens$c==chr,c('g','sn')]
  points(dens$sn, dens$g, type='l', col = c('red','dark red')[(chr%%2)+1], lwd=lwid)
  abline(h=chrborders[chr], col='dark grey',lty=2 , lwd=lwid)
}

axis(2, tckpos, seq(1,22,1), tick=FALSE, las=1)
axis(2, chrborders, rep('',23))
axis(3, c(0,20000), c(0,'20k'))
mtext('SNPs', side=2, line=13, adj=0.55, cex=10, outer=T)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')


#main plot
plot(c(chrborders[1], chrborders[23]), c(chrborders[1], chrborders[23]), xlab='', ylab='', xaxt='n', yaxt='n',col='white',xaxs = 'i', yaxs='i')
for(m in 1:(length(chrborders)-1)){
  i=ifelse(m%%2==0,1,2)
  for(s in 1:(length(chrborders)-1)){
    rect(chrborders[m], chrborders[s],chrborders[m+1], chrborders[s+1], col=color[i], border=NA)
    i=i+1
  }
}
print(date())
points(sigtab.cis[,1], sigtab.cis[,2], pch=20, col='green', cex=4)
points(sigtab.longrange[,1], sigtab.longrange[,2], pch=20, col='blue', cex=4)
points(sigtab.trans[,1], sigtab.trans[,2], pch=20, col='red', cex=4)
#points(sigtab.cis[,1], sigtab.cis[,2], pch=20, col='#006400', cex=2)
#points(sigtab.longrange[,1], sigtab.longrange[,2], pch=20, col='#FFA500', cex=2)
#points(sigtab.trans[,1], sigtab.trans[,2], pch=20, col='#1874CD', cex=2)
print(date())
box(lwd = 1.5)

#CpG density
plot(c(chrborders[1],chrborders[23]), c(0,max(res.dens$mn)), col='white', xlab='', ylab='', xaxt='n', yaxt='n',xaxs = 'i')
for(chr in 1:22){
  dens=res.dens[res.dens$c==chr,c('g','mn')]
  points(dens$g, dens$mn, type='l', col = c('red','dark red')[(chr%%2)+1], lwd=lwid)
  abline(v=chrborders[chr], col='dark grey',lty=2 , lwd=lwid)
}
axis(1, tckpos, seq(1,22,1), tick=FALSE, padj=1)
axis(1, chrborders, rep('',23))
axis(2, c(0,3000),c(0,'3k'), las=1)
mtext('CpGs', side=1, line= 23, adj=0.5, cex=10)

dev.off()


if(FALSE){
  load('/project/lolipop_b/METHQTL2/REPL2/cosmopairs_chessboard_070617.RData')
  test=abs(sigtab.longrange$gpos.cpg-sigtab.longrange$gpos.snp)
  png('/home/blehne/METHQTL/PLOTS/longrange_hist_070617.png', width = 600, height=600, units ='px')
  par(mar=c(6.5,6.5,5,3), mgp=c(4.2, 0.8, 0))
  hist(test, breaks=100, col='blue', main='Long-range Cis', xlab= 'Distance (mb)', xaxt='n', cex.lab=2.5, cex.main=3, cex.axis=1.5)
  axis(1, seq(0,250000000,50000000),seq(0,250,50), cex.axis=1.5)
  box()
  dev.off()
  
  png('/home/blehne/METHQTL/PLOTS/longrange_hist1_070617.png', width = 600, height=600, units ='px')
  par(mar=c(6.5,6.5,5,3), mgp=c(4.2, 0.8, 0))
  hist(test, breaks=250, col='blue', main='Long-range Cis', xlab= 'Distance (mb)', xaxt='n', xlim = c(0,15000000), cex.lab=2.5, cex.main=3, cex.axis=1.5)
  axis(1, seq(0,250000000,2500000),seq(0,250,2.5), cex.axis=1.5)
  box()
  dev.off()
  
  png('/home/blehne/METHQTL/PLOTS/longrange_hist3_070617.png', width = 900, height=600, units ='px')
  par(mar=c(6.5,6.5,5,3), mgp=c(4.2, 0.8, 0))
  r=hist(test, breaks=250, plot=F)
  barplot(log10(r$counts+1)[1:20],names.arg=1:20, col='blue', main='', xlab= 'Distance (mb)', cex.lab=2.5, cex.main=3, cex.axis=1.5, ylab='log10(Frequency)', ylim= c(0,8))
  box()
  dev.off()
  
  test=c(abs(sigtab.longrange$gpos.cpg-sigtab.longrange$gpos.snp), abs(sigtab.cis$gpos.cpg-sigtab.cis$gpos.snp))
  png('/home/blehne/METHQTL/PLOTS/longrange_hist4_120617.png', width = 900, height=600, units ='px')
  par(mar=c(6.5,6.5,5,3), mgp=c(4.2, 0.8, 0))
  r=hist(test, breaks=seq(0,250000000,100000), plot=F)
  barplot(log10(r$counts+1)[1:200],names.arg='', col='blue', main='', xlab= 'Distance (mb)', cex.lab=2.5, cex.main=3, cex.axis=1.5, ylab='log10(Frequency)', ylim= c(0,8))
  box()
  b=barplot(log10(r$counts+1)[1:300], plot=F)
  axis(1, b[seq(1,225,25)], seq(0,20,2.5))
  dev.off()
  
  # compare to trans
  load('/project/lolipop_b/METHQTL2/REPL2/cosmopairs_combined_151216.RData')
  sigtab.trans=cosmo[cosmo$cpg.chr!=cosmo$snp.chr,]
  res.trans=matrix(nrow=22, ncol=22)
  for(cpg.chr in 1:22){
    for(snp.chr in 1:22){
      res.trans[snp.chr, cpg.chr]=length(unique(sigtab.trans[sigtab.trans$cpg.chr==cpg.chr & sigtab.trans$snp.chr==snp.chr,'pair']))
    }
  }
  sigtab.cis=cosmo[cosmo$cpg.chr==cosmo$snp.chr,]
  thresholds=c(1,seq(1000000,20000000,1000000))
  res.cis=matrix(nrow=22, ncol=length(thresholds))
  cnames=NA
  for(chr in 1:22){
    for(c in 1:length(thresholds)){
      #res.cis=c(res.cis, length(unique(sigtab.cis[sigtab.cis$cpg.chr==chr,'pair'])))
      res.cis[chr, c]=length(unique(sigtab.cis[sigtab.cis$cpg.chr==chr & abs(sigtab.cis$snp.pos-sigtab.cis$cpg.pos)>thresholds[c],'pair']))
    }
  }
  save(res.cis, res.trans, file ='/project/lolipop_b/METHQTL2/REPL2/cisvstrans_140617.RData')
  
  
  # visualise all trans results 
  tmp1=apply(res.trans,2,median)
  names(tmp1)=seq(1,22,1)
  boxes1=res.trans[,as.numeric(as.character(names(sort(tmp1, decreasing=T))))]
  tmp2=apply(res.trans,1,median)
  names(tmp2)=seq(1,22,1)
  boxes2=res.trans[as.numeric(as.character(names(sort(tmp2, decreasing=T)))),]
  
  png('/home/blehne/METHQTL/PLOTS/longrange_boxplot_140617.png', width = 2500, height=1200, units ='px')
  par(mar=c(9,9,3,3), mgp=c(5.5,1.3,0))
  #boxplot(cbind(res.cis, boxes1, boxes2), col=c(rep('red', length(thresholds)), rep('dark green', 22), rep('dark blue', 22)), ylim = c(1,1000000), log='y', xaxt='n')
  boxplot(cbind(res.cis[,1:11], boxes1), col=c(rep('red', 11), rep('dark green', 22)), ylim = c(10,1000000), log='y', xaxt='n', xlab = 'Cis/Trans', ylab='Frequency', cex.lab=3, cex.axis=2)
  axis(1, seq(1,33,1), c('all cis', paste('>',seq(1,10,1), 'Mb',sep=''), paste('chr',as.numeric(as.character(names(sort(tmp1, decreasing=T)))), sep='')), cex.axis=1.4)
  dev.off()
  
  
  
  
  
  
  
  
  
}
