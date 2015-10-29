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
library(gplots)

## load custom functions
# source('code/custom_ggplot_settings.R')
# source('code/custom_heatmap.3.R')

## load OTU table and map file from setup
load('results/otu_setup/ger_otu_setup.RData')

## heatmap
## write dataframe with shorter taxa names
df.taxa <- data.frame(full = colnames(otu.g.prop),
                      short = '')

## subset to top 50 taxa
pick <- order(colSums(otu.g.prop), decreasing=TRUE)[1:50]
heatmap(otu.g.prop[, pick], trace = NULL)

####################################
## compute Bray-Curtis and Canberra distances on rarefied data
bc.g <- vegdist(otu.g.rare)
can.g <- vegdist(otu.g.rare, 'canberra')

####################################
## PCoA
pcoa.bc.g <- cmdscale(t(bc.g), eig = TRUE)
pcoa.can.g <- cmdscale(t(can.g), eig = TRUE)

####################################
## plot Bray-Curtis pcoa with ggplot
df.pcoa.bc.g <- as.data.frame(pcoa.bc.g$points)
colnames(df.pcoa.bc.g) <- c('PCoA1', 'PCoA2')
df.pcoa.bc.g$Description <- rownames(df.pcoa.bc.g)
df.pcoa.bc.g.all <- merge(df.pcoa.bc.g, ger.map)

## specify 12 colors for space types
mycol.12 <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', 
              '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#EEEE00', '#B15928')

## all samples by space type
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.bc.g + geom_point(size = 2) +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 2, vjust = 0) + 
  scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_sampletype.png')

## crack area NSF
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackAreaNSF))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_crackarea.png')

## chem = TCSavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_TCSavg.png')

## chem = TCCavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_TCCavg.png')

## chem = MePBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_MePBavg.png')

## chem = EtPBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_EtPBavg.png')

## chem = PrPBavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_PrPBavg.png')

## chem = BuBPavg  
gg.pcoa.bc.g <- ggplot(df.pcoa.bc.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.bc.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_bc_g_BuPBavg.png')

#################################
## plot Canberra pcoa with ggplot
df.pcoa.can.g <- as.data.frame(pcoa.can.g$points)
colnames(df.pcoa.can.g) <- c('PCoA1', 'PCoA2')
df.pcoa.can.g$Description <- rownames(df.pcoa.can.g)
df.pcoa.can.g.all <- merge(df.pcoa.can.g, ger.map)

## all samples by space type
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE))
gg.pcoa.can.g + geom_point(size = 2) +
  geom_text(aes(y = PCoA2 + 0.01, label = Description), size = 2, vjust = 0) + 
  scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_sampletype.png')

## crack area NSF
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = CrackAreaNSF))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_crackarea.png')

## chem = TCSavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCSavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_TCSavg.png')

## chem = TCCavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = TCCavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_TCCavg.png')

## chem = MePBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = MePBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_MePBavg.png')

## chem = EtPBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = EtPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_EtPBavg.png')

## chem = PrPBavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = PrPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_PrPBavg.png')

## chem = BuBPavg  
gg.pcoa.can.g <- ggplot(df.pcoa.can.g.all, aes(x = PCoA1, y = PCoA2, color = SpaceTypeBioBE, size = BuPBavg))
gg.pcoa.can.g + geom_point() + scale_color_manual(values = mycol.12)
ggsave('figures/pcoa_can_g_BuPBavg.png')
