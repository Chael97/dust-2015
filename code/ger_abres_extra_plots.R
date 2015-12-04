## Perform exploratory ecological analysis on taxonomic data
## 2015-12-03
## Roxana Hickey <roxana.hickey@gmail.com>

setwd('~/Documents/projects/gerlinger/')
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(labdsv)
library(plyr)
library(phyloseq)
library(gplots)

## load previously processed metagenomic data
load('results/otu_setup/ger_rm_contaminants_meta.RData')
load('~/Documents/projects/gerlinger/results/ger_shortbred_abres_plots.RData')

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
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control'))
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
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control'))
ggsave('figures/pcoa_bc_ardb_crackarea.png', width = 8, height = 6.5, units = 'in')

## chem = TCSavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Triclosan')
ggsave('figures/pcoa_bc_ardb_TCSavg.png', width = 8, height = 6.5, units = 'in')

## chem = TCCavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Trichlorocarbanilide')
ggsave('figures/pcoa_bc_ardb_TCCavg.png', width = 8, height = 6.5, units = 'in')

## chem = MePBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Methylparaben')
ggsave('figures/pcoa_bc_ardb_MePBavg.png', width = 8, height = 6.5, units = 'in')

## chem = EtPBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Ethylparaben')
ggsave('figures/pcoa_bc_ardb_EtPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = PrPBavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Propylparaben')
ggsave('figures/pcoa_bc_ardb_PrPBavg.png', width = 8, height = 6.5, units = 'in')

## chem = BuBPavg  
gg.pcoa.bc.ardb <- ggplot(df.pcoa.bc.ardb.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.ardb + 
  geom_point(data = df.pcoa.bc.ardb.all[df.pcoa.bc.ardb.all$ChemAvail == FALSE,], 
             aes(x = PCoA1, y = PCoA2), 
             size = 2, pch = 1) +
  geom_point() + 
  scale_color_manual(values = mycol.10, name = 'Space Type',
                     breaks = c('building.support', 'circulation', 'classroom', 'gym', 'laundry', 
                                'lockers', 'office', 'pool', 'restroom', 'neg.control')) +
  scale_size_continuous(range = c(3,8), name = 'Butylparaben')
ggsave('figures/pcoa_bc_ardb_BuPBavg.png', width = 8, height = 6.5, units = 'in')