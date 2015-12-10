## Perform exploratory ecological analysis on taxonomic data
## Last updated 2015-12-09
## Roxana Hickey <roxana.hickey@gmail.com>

setwd('~/Documents/projects/gerlinger/')
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(labdsv)
library(plyr)
library(phyloseq)
library(gplots)

## load OTU table and map file from post-contaminant filtering
# load('results/otu_setup/ger_otu_setup.RData')
load('results/otu_setup/ger_rm_contaminants.RData')

## Helpful reminders: 
## ger.nc.rare is the rarefied dataset (depth = 16,000) with contaminants and controls removed. samples are in rows.
## taxo.nc and consensus.nc are the taxonomic labels and short names for the rarefied dataset (ger.nc.rare)
## bc.nc is the Bray-curtis dissimilarity matrix based on ger.nc.rare
## can.nc is the Canberra dissimilarity matrix based on ger.nc.rare
## pcoa.bc.nc is the corresponding PCoA object to bc.nc

####################################
## heatmap
## make proportion table
ger.nc.prop <- prop.table(ger.nc, margin = 1)

## subset to top 25 taxa
# pick <- order(apply(ger.nc.prop, 2, median), decreasing=TRUE)[1:25]
pick <- order(colSums(ger.nc.prop), decreasing=TRUE)[1:25]

identical(row.names(ger.nc.prop), row.names(ger.map))
mycol <- mycol.9[match(ger.map$SpaceTypeBioBE, names(mycol.9))]

png('figures/ger_heatmap_top25_genus.png', width = 7, height = 6, res = 300, units = 'in', pointsize = 8)
par(xpd = TRUE)
heatmap.2(t(ger.nc.prop[,pick]), trace = 'none',
          col = colorRampPalette(brewer.pal(9, 'YlGnBu'))(100), margin = c(6,10),
          ColSideColors = mycol,
          density.info = 'none', 
          keysize = 1, key.title = '', key.xlab = 'Proportion',
          labCol = ger.map[row.names(ger.nc.prop), 'Description'],
          labRow = consensus.nc[colnames(ger.nc.prop[,pick])],
          xlab = 'Sample', ylab = 'Taxon')
legend(0.9, 1.1, bty = 'n', cex = 0.7, title = 'Space Type',
       legend = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                  'lockers', 'office', 'pool', 'restroom'),
       fill = mycol.9, border = 'white')
dev.off()

# png('figures/ger_heatmap_top25_genus_v2.png', width = 7, height = 6, res = 300, units = 'in', pointsize = 8)
# heatmap.2(t(ger.nc.prop[,pick]), trace = 'none',
#           col = colorRampPalette(brewer.pal(9, 'YlGnBu'))(100), margin = c(6,10),
#           # ColSideColors = mycol,
#           density.info = 'none', 
#           keysize = 1, key.title = '', key.xlab = 'Proportion',
#           labCol = ger.map[row.names(ger.nc.prop), 'Description'],
#           labRow = consensus.nc[colnames(ger.nc.prop[,pick])],
#           xlab = 'Sample', ylab = 'Taxon')
# dev.off()

####################################
## plot Bray-Curtis PCoA
df.pcoa.bc.nc <- as.data.frame(pcoa.bc.nc$points)
colnames(df.pcoa.bc.nc) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.nc$SampleID <- rownames(df.pcoa.bc.nc)
df.pcoa.bc.nc.all <- merge(df.pcoa.bc.nc, ger.map)

## set theme_bw
theme_set(theme_bw())

## all samples by space type
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.nc + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_bc_g_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Bray-Curtis NMDS
df.nmds.bc.nc <- as.data.frame(nmds.bc.nc$points)
colnames(df.nmds.bc.nc) <- c('NMDS1', 'NMDS2')
df.nmds.bc.nc$SampleID <- rownames(df.nmds.bc.nc)
df.nmds.bc.nc.all <- merge(df.nmds.bc.nc, ger.map)

## all samples by space type
gg.nmds.bc.nc <- ggplot(df.nmds.bc.nc.all, aes(x = NMDS1, y = NMDS2, color = SpaceTypeBioBE))
gg.nmds.bc.nc + geom_point(size = 3) +
  # ggtitle('NMDS on Bray-Curtis') +
  geom_text(aes(y = NMDS2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/nmds_bc_g_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Bray-Curtis pcoa with crack/chem data
## crack area NSF
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackAreaNSF))
gg.pcoa.bc.nc +   
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$CrackNSFAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1, na.rm = FALSE) +
  geom_point() + 
  scale_size_continuous(range = c(3,8), name = 'Crack Area') +
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_bc_g_crackarea.png', width = 8, height = 6.5, units = 'in')

## chem = TCSavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclosan')
ggsave('figures/pcoa_bc_g_TCSavg.png', width = 8, height = 6.5, units = 'in')

## chem = TCCavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclocarban')
ggsave('figures/pcoa_bc_g_TCCavg.png', width = 8, height = 6.5, units = 'in')

## chem = MePBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Methylparaben')
ggsave('figures/pcoa_bc_g_MePBavg.png', width = 8, height = 6.5, units = 'in')

## chem = EtPBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Ethylparaben')
ggsave('figures/pcoa_bc_g_EtPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = PrPBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Propylparaben')
ggsave('figures/pcoa_bc_g_PrPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = BuBPavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Butylparaben')
ggsave('figures/pcoa_bc_g_BuPBavg.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Canberra PCoA
df.pcoa.can.nc <- as.data.frame(pcoa.can.nc$points)
colnames(df.pcoa.can.nc) <- c('PCoA1', 'PCoA2')
df.pcoa.can.nc$SampleID <- rownames(df.pcoa.can.nc)
df.pcoa.can.nc.all <- merge(df.pcoa.can.nc, ger.map)

## set theme_bw
theme_set(theme_bw())

## all samples by space type
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.can.nc + geom_point(size = 3) +
  # ggtitle('PCoA on Canberra') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_can_g_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Canberra NMDS
df.nmds.can.nc <- as.data.frame(nmds.can.nc$points)
colnames(df.nmds.can.nc) <- c('NMDS1', 'NMDS2')
df.nmds.can.nc$SampleID <- rownames(df.nmds.can.nc)
df.nmds.can.nc.all <- merge(df.nmds.can.nc, ger.map)

## all samples by space type
gg.nmds.can.nc <- ggplot(df.nmds.can.nc.all, aes(x = NMDS1, y = NMDS2, color = SpaceTypeBioBE))
gg.nmds.can.nc + geom_point(size = 3) +
  # ggtitle('NMDS on Canberra') +
  geom_text(aes(y = NMDS2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/nmds_can_g_sampletype.png', width = 8, height = 6.5, units = 'in')

####################################
## plot Canberra pcoa with crack/chem data
## crack area NSF
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackAreaNSF))
gg.pcoa.can.nc +   
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$CrackNSFAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1, na.rm = FALSE) +
  geom_point() + 
  scale_size_continuous(range = c(3,8), name = 'Crack Area') +
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom'))
ggsave('figures/pcoa_can_g_crackarea.png', width = 8, height = 6.5, units = 'in')

## chem = TCSavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclosan')
ggsave('figures/pcoa_can_g_TCSavg.png', width = 8, height = 6.5, units = 'in')

## chem = TCCavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Triclocarban')
ggsave('figures/pcoa_can_g_TCCavg.png', width = 8, height = 6.5, units = 'in')

## chem = MePBavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Methylparaben')
ggsave('figures/pcoa_can_g_MePBavg.png', width = 8, height = 6.5, units = 'in')

## chem = EtPBavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Ethylparaben')
ggsave('figures/pcoa_can_g_EtPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = PrPBavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Propylparaben')
ggsave('figures/pcoa_can_g_PrPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = BuBPavg  
gg.pcoa.can.nc <- ggplot(df.pcoa.can.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.can.nc + 
  geom_point(data = df.pcoa.can.nc.all[df.pcoa.can.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.9, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom')) +
  scale_size_continuous(range = c(3,8), name = 'Butylparaben')
ggsave('figures/pcoa_can_g_BuPBavg.png', width = 8, height = 6.5, units = 'in')