library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scales)
library(cowplot)
theme_set(theme_cowplot() + background_grid(major="x"))

setwd("C:/Users/Johann Hawe/Work/data_transfer/hmgu/meqtl_paper")
load("data/cosmopairs_combined_151216.RData")

cosmo_trans <- cosmo %>% mutate(is_trans = snp.chr != cpg.chr) %>% 
  filter(is_trans)
cosmo_trans <- 
  cosmo_trans %>% 
  group_by(snp, snp.chr, snp.pos) %>% 
  summarise(total_associations = n()) %>% 
  ungroup() %>% 
  arrange(desc(total_associations))

# get the hg19 chromosome definitions
hg19info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

cosmo_ranges <- with(cosmo_trans,
                     GRanges(
                       paste0("chr", snp.chr),
                       IRanges(snp.pos, width = 2),
                       name = snp,
                       nassocs = total_associations
                     ))

# we need to convert chromsome positions to proper xaxis positions

# get chromosome information (lengths) and create bins
chrs = paste0("chr", 1:22)
chrlen <- seqlengths(hg19info)
chrlen <- chrlen[chrs]
genome_bins <- tileGenome(chrlen,
                          tilewidth = 1e7,
                          cut.last.tile.in.chrom = T)

# define the breaks for plotting boundaries
breaks <- table(seqnames(genome_bins))
for (i in 2:length(breaks)) {
  breaks[i] <- breaks[i - 1] + breaks[i]
}

# get bin overlaps of our cosmo ranges
x_bin <- subjectHits(findOverlaps(cosmo_ranges, genome_bins))
top_percent_associations <- quantile(cosmo_trans$total_associations, 0.99)

# finalize plot data
plot_data <- bind_cols(cosmo_trans, xbin = x_bin) %>%
  mutate(is_top = total_associations > top_percent_associations) %>%
  distinct(xbin, total_associations, is_top)

# helper to squish y axis
squish_trans <- function(from, to, factor) { 
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- which(x > from & x < to)
    ito <- which(x >= to)
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- which(x > from & x < from + (to - from)/factor)
    ito <- which(x >= from + (to - from)/factor)
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squish_and_reverse",
                   transform = trans, 
                   inverse = inv))
}

gp <- ggplot(plot_data,
       aes(y = total_associations, x = xbin, color = is_top)) +
  geom_point() + 
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(trans = squish_trans(60, 280, 20),
                     labels = c(0,20,40,60,280,300),
                     breaks = c(0,20,40,60,280,300)) +
  scale_x_continuous(breaks = c(0, as.vector(breaks)),
                     labels = c(names(breaks), "")) +
  
  labs(x = "Chromosome", y = "Associated CpGs per SNP") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ), legend.position = "none")

save_plot("manhattan.pdf", gp, base_asp = 2, ncol=1.5, nrow=1.5)
