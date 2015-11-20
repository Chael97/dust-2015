## Perform exploratory ecological analysis on taxonomic data
## 2015-10-27
## Roxana Hickey <roxana.hickey@gmail.com>

setwd('~/Documents/projects/gerlinger/')
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(labdsv)
library(plyr)
library(phyloseq)
library(gplots)

## load custom functions
# source('code/custom_ggplot_settings.R')
# source('code/custom_heatmap.3.R')

## load OTU table and map file from post-contaminant filtering
# load('results/otu_setup/ger_otu_setup.RData')
load('results/otu_setup/ger_rm_contaminants.RData')

## Helpful reminders: 
## ger.nc.6500 is the rarefied dataset with contaminants removed, controls retained. samples are in rows.
## taxo.nc and consensus.nc are the taxonomic labels and short names for the rarefied dataset (ger.nc.6500)
## bc.nc is the Bray-curtis dissimilarity matrix based on ger.nc.6500
## pcoa.bc.nc is the corresponding PCoA object to bc.nc

####################################
## heatmap
## make proportion table
ger.nc.prop <- prop.table(ger.nc, margin = 1)

## subset to top 25 taxa
# pick <- order(apply(ger.nc.prop, 2, median), decreasing=TRUE)[1:25]
pick <- order(colSums(ger.nc.prop), decreasing=TRUE)[1:25]

png('figures/ger_heatmap_top25_genus.png', width = 6, height = 6, res = 300, units = 'in', pointsize = 8)
heatmap.2(ger.nc.prop[,pick], trace = 'none',
          col = colorRampPalette(brewer.pal(9, 'YlGnBu'))(100), margin = c(10,10), 
          keysize = 1, density.info = 'none', key.title = '', key.xlab = 'Proportion',
          labRow = ger.map[row.names(ger.nc.prop), 'SpaceType'],
          labCol = consensus.nc[colnames(ger.nc.prop[,pick])])
dev.off()

####################################
## plot Bray-Curtis PCoA
df.pcoa.bc.nc <- as.data.frame(pcoa.bc.nc$points)
colnames(df.pcoa.bc.nc) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.nc$SampleID <- rownames(df.pcoa.bc.nc)
df.pcoa.bc.nc.all <- merge(df.pcoa.bc.nc, ger.map)

## specify 12 colors for space types
mycol.12 <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', 
              '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#EEEE00', '#B15928')

## set theme_bw
theme_set(theme_bw())

## all samples by space type
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.nc + geom_point(size = 3) +
  # ggtitle('PCoA on Bray-Curtis') +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 3, vjust = 0) + 
  scale_color_manual(values = mycol.12)
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
  scale_color_manual(values = mycol.12)
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
  scale_size_continuous(range = c(3,8)) +
  scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_crackarea.png', width = 8, height = 6.5, units = 'in')

## chem = TCSavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_TCSavg.png', width = 8, height = 6.5, units = 'in')

## chem = TCCavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_TCCavg.png', width = 8, height = 6.5, units = 'in')

## chem = MePBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_MePBavg.png', width = 8, height = 6.5, units = 'in')

## chem = EtPBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_EtPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = PrPBavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_PrPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = BuBPavg  
gg.pcoa.bc.nc <- ggplot(df.pcoa.bc.nc.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.nc + 
  geom_point(data = df.pcoa.bc.nc.all[df.pcoa.bc.nc.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.12) + 
  scale_size_continuous(range = c(3,8))
ggsave('figures/pcoa_bc_g_BuPBavg.png', width = 8, height = 6.5, units = 'in')
