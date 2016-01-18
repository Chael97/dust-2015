## Perform exploratory analysis on antibiotic resistome data
## Roxana Hickey <roxana.hickey@gmail.com>
## Last updated 2015-12-09

setwd('~/Documents/projects/dust_2015/')
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(labdsv)
library(plyr)
library(phyloseq)
library(gplots)

## load previously processed metagenomic data
load('results/otu_setup/ger_rm_contaminants_meta.RData')
load('~/Documents/projects/dust_2015/results/ger_shortbred_abres_plots.RData')

## pcoa based on ARDB genes
bc.ardb <- vegdist(t(ger.ardb))
pcoa.bc.ardb <- cmdscale(bc.ardb, eig = TRUE)

plot(pcoa.bc.ardb$points, type = 'p')

## plot Bray-Curtis PCoA
df.pcoa.bc.ardb <- as.data.frame(pcoa.bc.ardb$points)
colnames(df.pcoa.bc.ardb) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.ardb$SampleID <- rownames(df.pcoa.bc.ardb)
df.pcoa.bc.ardb.all <- merge(df.pcoa.bc.ardb, ger.meta.map)

eigenvals(pcoa.bc.ardb)/sum(eigenvals(pcoa.bc.ardb))
r2.pcoa.bc.ardb.1 <- round(eigenvals(pcoa.bc.ardb)[1]/sum(eigenvals(pcoa.bc.ardb)) * 100, 1)
r2.pcoa.bc.ardb.2 <- round(eigenvals(pcoa.bc.ardb)[2]/sum(eigenvals(pcoa.bc.ardb)) * 100, 1)

## set theme_bw
theme_set(theme_bw(base_size = 14))

## all samples by space type
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.ardb + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  xlab(paste0('PCoA 1 (', r2.pcoa.bc.ardb.1, '%)')) +
  ylab(paste0('PCoA 2 (', r2.pcoa.bc.ardb.2, '%)')) +
  theme(panel.grid.minor = element_blank())
ggsave('figures/pcoa_bc_ardb_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## pcoa based on CARD genes
bc.card <- vegdist(t(ger.card))
pcoa.bc.card <- cmdscale(bc.card, eig = TRUE)

plot(pcoa.bc.card$points, type = 'p')

eigenvals(pcoa.bc.card)/sum(eigenvals(pcoa.bc.card))
r2.pcoa.bc.card.1 <- round(eigenvals(pcoa.bc.card)[1]/sum(eigenvals(pcoa.bc.card)) * 100, 1)
r2.pcoa.bc.card.2 <- round(eigenvals(pcoa.bc.card)[2]/sum(eigenvals(pcoa.bc.card)) * 100, 1)

## plot Bray-Curtis PCoA
df.pcoa.bc.card <- as.data.frame(pcoa.bc.card$points)
colnames(df.pcoa.bc.card) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.card$SampleID <- rownames(df.pcoa.bc.card)
df.pcoa.bc.card.all <- merge(df.pcoa.bc.card, ger.meta.map)

## all samples by space type
gg.pcoa.bc.card <- ggplot(df.pcoa.bc.card.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.card + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  xlab(paste0('PCoA 1 (', r2.pcoa.bc.ardb.1, '%)')) +
  ylab(paste0('PCoA 2 (', r2.pcoa.bc.ardb.2, '%)')) +
  theme(panel.grid.minor = element_blank())
ggsave('figures/pcoa_bc_card_sampletype.png', width = 8, height = 6.5, units = 'in')

save.image('results/ger_abres_extra_plots.RData')