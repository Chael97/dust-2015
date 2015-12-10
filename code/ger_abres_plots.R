## Summarize Gerlinger metagenomics ShortBred data
## Roxana Hickey <roxana.hickey@gmail.com
## Last updated 2015-12-09

setwd('~/Documents/projects/gerlinger/')
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(MASS)
library(scales)

## load OTU table and map file from post-contaminant filtering
load('results/otu_setup/ger_rm_contaminants_meta.RData')

## Helpful reminders: 
## ger.meta.nc is the MetaPhlAn dataset with Propis and control removed. samples are in rows.
## taxo.nc and consensus.nc are the taxonomic labels and short names for the dataset (ger.meta.nc)
## bc.meta.nc is the Bray-curtis dissimilarity matrix based on ger.meta.nc
## pcoa.bc.meta.nc is the corresponding PCoA object to bc.meta.nc

####################################
## combined MBTA/human ARDB data from Huttenhower group
hut.ardb <- read.table('data/shortbred_ardb/ARDB_Data_comparison.txt',
                       sep = '\t', header = T, row.names = 1)

## Huttenhower metadata
hut.cat <- read.table('data/shortbred_ardb/categories_metadata_v2_151130.txt',
                      sep = '\t', header = 1, row.names = NULL)

## Gerlinger ARDB data
ger.ardb <- read.table('data/shortbred_ardb/151203_gerlinger_shortbred_merged_ardb_barcode_sampleid.txt',
                       sep = '\t', header = T, row.names = 1)

## remove zero-sum samples (no sequence data)
ger.ardb <- ger.ardb[,colSums(ger.ardb) > 0]

## check if all samples in master sample ID key
colnames(ger.ardb) %in% sid.key$metaphlan_barcode
colnames(ger.ardb)[!(colnames(ger.ardb) %in% sid.key$metaphlan_barcode)]

## remove 'undetermined' sample
ger.ardb[,'Undetermined'] <- NULL
colnames(ger.ardb) %in% sid.key$metaphlan_barcode

## replace sample IDs to match 16S data
colnames(ger.ardb) <- sid.key$uparse_sampleid[match(colnames(ger.ardb), sid.key$metaphlan_barcode)]

## remove duplicate samples
ger.ardb <- ger.ardb[, colnames(ger.ardb) %in% rownames(ger.meta)]

## read in key with more informative AbRes class labels
ardb.fam.id <- as.data.frame(read.table('data/shortbred_ardb/ARDB_family_descriptions.txt', 
                                        header = T, sep = '\t'))

## Gerlinger metadata
ger.cat <- data.frame('SampleID' = colnames(ger.ardb),
                      'SampleType' = rep('Indoor', length(colnames(ger.ardb))),
                      'Env' = rep('Gerlinger', length(colnames(ger.ardb))))

## check that datasets match
identical(rownames(ger.ardb), rownames(hut.ardb)) # TRUE

## combine datasets and sort
all.ardb <- cbind(ger.ardb, hut.ardb)
all.cat <- rbind(ger.cat, hut.cat)

colnames(all.ardb) %in% all.cat$SampleID

## simplify all antibiotics to one summed RPKM measurement
all.ardb.sum <- colSums(all.ardb)
all.cat$AbRes_RPKM <- all.ardb.sum[match(names(all.ardb.sum), all.cat$SampleID)]

## set theme_bw
theme_set(theme_bw())

## violin plot of 3 studies
gg.ardb.vio <- ggplot(all.cat, aes(x = Env, y = AbRes_RPKM, color = Env, fill = Env))
gg.ardb.vio + geom_violin(adjust = 0.75, alpha = 0.3) +
  geom_boxplot(width = 0.1, aes(fill = Env), color = 'black', outlier.color = NA) +
  stat_summary(fun.y = median, geom = 'point', fill = 'white', shape = 21, size = 2.5) +
  scale_color_manual(values = c('darkmagenta', '#0049FF', '#00A4DE'),
                     limits = c('Gerlinger',                                                               
                                'Built_Environment',
                                'Human')) +
  scale_fill_manual(values = c('darkmagenta', '#0049FF', '#00A4DE'),
                    limits = c('Gerlinger',                                                               
                               'Built_Environment',
                               'Human')) +
  scale_x_discrete(limits = c('Gerlinger',                                                               
                              'Built_Environment',
                              'Human'), 
                   labels = c('Present Study',                                                               
                              'Other Built\nEnvironment Studies',
                              'Human Microbiome\nStudies')) +
  scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x))) +
  ylab('Total antibiotic resistance gene families (RPKM)') + 
  theme(legend.position = 'none', axis.title.x = element_blank())
ggsave('figures/ger_mbta_hmp_ardb_violin_log.png', width = 8, height = 8)

#####################################################
## NOTE TO SELF: Update AbResFam with informative IDs
#####################################################

## prep data in long format
all.ardb.tp <- as.data.frame(t(all.ardb))
all.ardb.tp$SampleID <- rownames(all.ardb.tp)
all.ardb.tp$Env <- factor(all.cat$Env[match(all.cat$SampleID, all.ardb.tp$SampleID)], 
                          levels = c('Human', 'Built_Environment', 'Gerlinger'))

all.ardb.lg <- melt(all.ardb.tp, id.vars = c('SampleID', 'Env'))
colnames(all.ardb.lg) <- c('SampleID', 'Env', 'AbResFam', 'RPKM')

## dot plot of all ab res gene families
gg.ardb.dot <- ggplot(all.ardb.lg, aes(x = AbResFam, y = RPKM, color = Env, 
                                       alpha = Env, shape = Env))
gg.ardb.dot + geom_point(position = 'jitter', aes(order = rev(seq(1, nrow(all.ardb.lg))))) +
  scale_color_manual(name = 'Study',
                     limits = c('Gerlinger', 'Built_Environment', 'Human'),
                     labels = c('Present Study',
                                'Other Built Environments',
                                'Human Microbiome'),
                     values = c('darkmagenta', '#0049FF', '#00A4DE')) +
  scale_shape_manual(name = 'Study',
                     limits = c('Gerlinger', 'Built_Environment', 'Human'),
                     labels = c('Present Study',
                                'Other Built Environments',
                                'Human Microbiome'),
                     values = c(17, 2, 1)) +
  scale_alpha_discrete(name = 'Study',
                       limits = c('Gerlinger', 'Built_Environment', 'Human'),
                       labels = c('Present Study',
                                  'Other Built Environments',
                                  'Human Microbiome'),
                       range = rev(c(0.1, 1))) +
  xlab('Antibiotic Resistance Gene Families') +
  scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        legend.position = 'bottom')
ggsave('figures/ger_mbta_hmp_ardb_dot_log.png', width = 8, height = 6)

## dot plot of all ab res gene families in Gerlinger found in 2+ samples
ger.ardb.1plus <- as.data.frame(t(ger.ardb[rowSums(ger.ardb) > 0,]))
ger.ardb.2plus  <- ger.ardb.1plus[,colSums(ger.ardb.1plus != 0) >= 2]
ger.ardb.2plus$SampleID <- rownames(ger.ardb.2plus)
ger.ardb.2plus.lg <- melt(ger.ardb.2plus, id.vars = 'SampleID')
colnames(ger.ardb.2plus.lg) <- c('SampleID', 'AbResFam', 'RPKM')

## add more informative AbRes names
ger.ardb.2plus.lg$AbResFamClass <- ardb.fam.id$Class[match(ger.ardb.2plus.lg$AbResFam, ardb.fam.id$Family)]
ger.ardb.2plus.lg$AbResFamClass2 <- paste0(ger.ardb.2plus.lg$AbResFamClass, ' (', ger.ardb.2plus.lg$AbResFam, ')')

ardb.fam.id.2plus <- ardb.fam.id[ardb.fam.id$Family %in% ger.ardb.2plus.lg$AbResFam,]

## plot wrt AbRes family
gg.ardb.ger.dot.f <- ggplot(ger.ardb.2plus.lg, aes(x = AbResFam, y = RPKM, color = AbResFamClass))
gg.ardb.ger.dot.f + geom_jitter(alpha = 0.7, position = position_jitter(width = 0.1),
                                aes(shape = AbResFamClass)) +
  xlab('ARDB Gene Family') +
  scale_x_discrete(limits = ardb.fam.id.2plus$Family[order(ardb.fam.id.2plus$Class)]) +
  scale_shape_manual(values = c(16,15,17,16,15,17,16,15,17,16,15,17,16)) +
  theme_bw() +
  guides(color = guide_legend(title = 'ARDB Gene Class'),
         shape = guide_legend(title = 'ARDB Gene Class')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('figures/ger_ardb_dot_family.png', width = 12, height = 7)

## plot wrt AbRes class
gg.ardb.ger.dot.c <- ggplot(ger.ardb.2plus.lg, aes(x = AbResFamClass, y = RPKM))
gg.ardb.ger.dot.c + geom_jitter(alpha = 0.5, color = 'darkmagenta', shape = 17, 
                                position = position_jitter(width = 0.2),
                                aes(order = rev(seq(1, nrow(ger.ardb.2plus.lg))))) +
  xlab('ARDB Gene Class') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('figures/ger_ardb_dot_class.png', width = 8, height = 6)

## summary stats for RPKM
summary(all.cat$AbRes_RPKM[all.cat$Env == 'Gerlinger'])
summary(all.cat$AbRes_RPKM[all.cat$Env == 'Built_Environment'])
summary(all.cat$AbRes_RPKM[all.cat$Env == 'Human'])

################################
## Gerlinger CARD data
ger.card <- read.table('data/shortbred_card/151207_gerlinger_shortbred_merged_card_barcode_sampleid.txt',
                       sep = '\t', header = T, row.names = 1)

## remove zero-sum samples (no sequence data)
ger.card <- ger.card[ , colSums(ger.card) > 0]

## check if all samples in master sample ID key
colnames(ger.card) %in% sid.key$metaphlan_barcode
colnames(ger.card)[!(colnames(ger.card) %in% sid.key$metaphlan_barcode)]

## remove 'undetermined' sample
ger.card[,'Undetermined'] <- NULL
colnames(ger.card) %in% sid.key$metaphlan_barcode

## replace sample IDs to match 16S data
colnames(ger.card) <- sid.key$uparse_sampleid[match(colnames(ger.card), sid.key$metaphlan_barcode)]

## remove duplicate samples
ger.card <- ger.card[, colnames(ger.card) %in% rownames(ger.meta)]

save.image('~/Documents/projects/gerlinger/results/ger_shortbred_abres_plots.RData')