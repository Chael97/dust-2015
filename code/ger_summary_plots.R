## Perform exploratory ecological analysis on taxonomic data
## 2015-10-27
## Roxana Hickey <roxana.hickey@gmail.com>

setwd('~/Documents/gerlinger/')
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(labdsv)
library(plyr)
library(phyloseq)

## load custom functions
source('code/custom_ggplot_settings.R')
source('code/custom_heatmap.3.R')

## load OTU table and map file from setup
load('results/otu_setup/ger_otu_setup.RData')

## heatmap
## subset to top 50 taxa
# pick <- order(colSums(otu.g), decreasing=TRUE)[1:25]
# heatmap(otu.g[,pick], trace = NULL, margins = c(15,5))

## bar plots

## compute Canberra and Bray-Curtis distances on rarefied data
can.g <- vegdist(otu.g.rare, 'canberra')
bc.g <- vegdist(otu.g.rare)

## PCoA
pcoa.can.g <- cmdscale(t(can.g), eig = TRUE)
# png('figures/pcoa_genus_can.png', width = 6, height = 6, units = 'in', res = 300)
# ordiplot(scores(pcoa.can.g)[,c(1,2)], type = 't', 
         # main = 'PCoA with genera (Canberra distance)',
         # xlab = 'PCoA 1', ylab = 'PCoA 2')
# dev.off()

pcoa.bc.g <- cmdscale(t(bc.g), eig = TRUE)
# png('figures/pcoa_genus_bc.png', width = 6, height = 6, units = 'in', res = 300)
# ordiplot(scores(pcoa.bc.g)[,c(1,2)], type = 't', 
         # main = 'PCoA with genera (Bray-Curtis dissimilarity)',
         # xlab = 'PCoA 1', ylab = 'PCoA 2')
# dev.off()

####################################
## plot Bray-Curtis pcoa with ggplot
df.pcoa.bc.g <- as.data.frame(pcoa.bc.g$points)
colnames(df.pcoa.bc.g) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.g$Description <- rownames(df.pcoa.bc.g)
df.pcoa.bc.g.all <- merge(df.pcoa.bc.g, ger.map)

## all samples by space type
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.g + geom_point(size = 3) + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_sampletype.png')

## crack area
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackArea))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_crackarea.png')

## chem = TCSavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_TCSavg.png')

## chem = TCCavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_TCCavg.png')

## chem = MePBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_MePBavg.png')

## chem = EtPBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_EtPBavg.png')

## chem = PrPBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_PrPBavg.png')

## chem = BuBPavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.bc.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_bc_g_BuPBavg.png')

#################################
## plot Canberra pcoa with ggplot
df.pcoa.can.g <- as.data.frame(pcoa.can.g$points)
colnames(df.pcoa.can.g) <- c('PCoA1', 'PCoA2')
df.pcoa.can.g$Description <- rownames(df.pcoa.can.g)
df.pcoa.can.g.all <- merge(df.pcoa.can.g, ger.map)

## all samples by space type
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.can.g + geom_point(size = 3) + 
  scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_sampletype.png')

## crack area  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackArea))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_crackarea.png')

## chem = TCSavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_TCSavg.png')

## chem = TCCavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_TCCavg.png')

## chem = MePBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_MePBavg.png')

## chem = EtPBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_EtPBavg.png')

## chem = PrPBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_PrPBavg.png')

## chem = BuBPavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = brewer.pal(length(unique(df.pcoa.can.g.all$SpaceTypeBioBE)), 'Paired'))
ggsave('figures/pcoa_can_g_BuPBavg.png')
