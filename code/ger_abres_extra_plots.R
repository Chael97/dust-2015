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

## set theme_bw
theme_set(theme_bw())

## all samples by space type
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.ardb + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_bc_ardb_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Bray-Curtis pcoa with crack/chem data
## crack area NSF
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackAreaNSF))
gg.pcoa.bc.ardb +   
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$CrackNSFAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1, na.rm = FALSE) +
  geom_point() + 
  scale_size_continuous(range = c(3,8), name = 'Crack Area') +
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_bc_ardb_crackarea.png', width = 8, height = 6.5, units = 'in')

## chem = TCSavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclosan')
ggsave('figures/pcoa_bc_ardb_TCSavg.png', width = 8, height = 6.5, units = 'in')

## chem = TCCavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclocarban')
ggsave('figures/pcoa_bc_ardb_TCCavg.png', width = 8, height = 6.5, units = 'in')

## chem = MePBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Methylparaben')
ggsave('figures/pcoa_bc_ardb_MePBavg.png', width = 8, height = 6.5, units = 'in')

## chem = EtPBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Ethylparaben')
ggsave('figures/pcoa_bc_ardb_EtPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = PrPBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Propylparaben')
ggsave('figures/pcoa_bc_ardb_PrPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = BuBPavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Butylparaben')
ggsave('figures/pcoa_bc_ardb_BuPBavg.png', width = 8, height = 6.5, units = 'in')

## plot RPKM ~ [chem]
ger.cat$AbRes_RPKM <- colSums(ger.ardb)
ger.cat.meta.map <- merge(ger.cat, ger.meta.map)
ger.ardb.chem <- ger.cat.meta.map[,c('SampleID', 'AbRes_RPKM', 'TCSavg', 'TCCavg', 
                                     'MePBavg', 'EtPBavg', 'PrPBavg', 'BuPBavg')]
ger.ardb.chem.lg <- melt(ger.ardb.chem, id.vars = c('SampleID', 'AbRes_RPKM'))
colnames(ger.ardb.chem.lg) <- c('SampleID', 'AbRes_RPKM', 'Chemical', 'AvgValue')

gg.ardb.chem.line <- ggplot(ger.ardb.chem.lg, aes(y = AbRes_RPKM, x = AvgValue, group = Chemical, color = Chemical))
gg.ardb.chem.line + geom_point() +
  stat_smooth(method = 'loess', aes(fill = Chemical)) +
  facet_wrap( ~ Chemical, nrow = 2, scales = 'free') +
  ylab('Total antibiotic resistance gene families (RPKM)') +
  xlab('Mean chemical concentration (ng/g)') +
  scale_color_manual(values = c('#0000A9', '#0049FF', '#00A4DE', 'turquoise', 'chartreuse3', '#FFD701')) +
  scale_fill_manual(values = c('#0000A9', '#0049FF', '#00A4DE', 'turquoise', 'chartreuse3', '#FFD701')) +
  theme(legend.position = 'none')
ggsave('figures/chem_vs_ardb_loess.png', width = 8, height = 5, units = 'in')

gg.ardb.chem.line + geom_point() +
  stat_smooth(method = 'lm', aes(fill = Chemical)) +
  facet_wrap( ~ Chemical, nrow = 2, scales = 'free') +
  ylab('Total antibiotic resistance gene families (RPKM)') +
  xlab('Mean chemical concentration (ng/g)') +
  scale_color_manual(values = c('#0000A9', '#0049FF', '#00A4DE', 'turquoise', 'chartreuse3', '#FFD701')) +
  scale_fill_manual(values = c('#0000A9', '#0049FF', '#00A4DE', 'turquoise', 'chartreuse3', '#FFD701')) +
  theme(legend.position = 'none')
ggsave('figures/chem_vs_ardb_linear.png', width = 8, height = 5, units = 'in')

####################################
## pcoa based on CARD genes
bc.card <- vegdist(t(ger.card))
pcoa.bc.card <- cmdscale(bc.card, eig = TRUE)

plot(pcoa.bc.card$points, type = 'p')

## plot Bray-Curtis PCoA
df.pcoa.bc.card <- as.data.frame(pcoa.bc.card$points)
colnames(df.pcoa.bc.card) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.card$SampleID <- rownames(df.pcoa.bc.card)
df.pcoa.bc.card.all <- merge(df.pcoa.bc.card, ger.meta.map)

## set theme_bw
theme_set(theme_bw())

## all samples by space type
gg.pcoa.bc.card <- ggplot(df.pcoa.bc.card.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.card + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_bc_card_sampletype.png', width = 8, height = 6.5, units = 'in')

save.image('results/ger_abres_extra_plots.RData')