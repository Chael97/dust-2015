## Comparison of QIIME, UPARSE and shotgun profiling of Gerlinger samples
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-09-25

## The goal of this analysis is to compare Erica's Gerlinger (dust) data processed by QIIME (16S rRNA), 
## UPARSE (16S rRNA) and MetaPhlAn (shotgun metagenomes) and determine whether they give different 
## profiles of community composition.

## initial setup
library(vegan)

## read in master sample ID key
sample.key <- read.table('data/sampleid_master_key.txt', header = T, sep = '\t')

## exclude the '5 µL' samples (poor sample quality)
sample.key.red <- sample.key[grep('05uL', sample.key$master_sampleid, invert = TRUE),]

## taxonomic levels key
tax.key <- data.frame(level = c('phylum', 'class', 'order', 'family', 'genus'),
                      label = c('L2' ,'L3', 'L4', 'L5', 'L6'))

## read in QIIME-16S taxon tables
files.qiime <- c('data/tax_16S_qiime/prop_tables/otu_table_merged_meta_L2.txt',
                 'data/tax_16S_qiime/prop_tables/otu_table_merged_meta_L3.txt',
                 'data/tax_16S_qiime/prop_tables/otu_table_merged_meta_L4.txt',
                 'data/tax_16S_qiime/prop_tables/otu_table_merged_meta_L5.txt',
                 'data/tax_16S_qiime/prop_tables/otu_table_merged_meta_L6.txt')

qiime.tax <- lapply(files.qiime, function(x) read.table(x, header=TRUE, sep='\t', row.names=1))

names(qiime.tax) <- c('phylum', 'class', 'order', 'family', 'genus')

## rename columns with master sample IDs
qiime.tax <- lapply(seq_along(qiime.tax), function(x){
  colnames(qiime.tax[[x]]) <- sample.key$master_sampleid[match(colnames(qiime.tax[[x]]),
                                                               sample.key$qiime_sampleid)]
  qiime.tax[[x]]})

names(qiime.tax) <- c('phylum', 'class', 'order', 'family', 'genus')

## reduce dataset to exclude '5 µL' samples
qiime.tax.red <- lapply(seq_along(qiime.tax), function(x){
  pick <- colnames(qiime.tax[[x]])[colnames(qiime.tax[[x]]) %in% sample.key.red$master_sampleid]
  qiime.tax[[x]][,pick]})

names(qiime.tax.red) <- c('phylum', 'class', 'order', 'family', 'genus')

## read in UPARSE-16S taxon tables
files.uparse <- c('data/tax_16S_uparse/prop_tables/fiveuL1201_demuxed_tax_normalized_otus_merged_meta_L2.txt',
                  'data/tax_16S_uparse/prop_tables/fiveuL1201_demuxed_tax_normalized_otus_merged_meta_L3.txt',
                  'data/tax_16S_uparse/prop_tables/fiveuL1201_demuxed_tax_normalized_otus_merged_meta_L4.txt',
                  'data/tax_16S_uparse/prop_tables/fiveuL1201_demuxed_tax_normalized_otus_merged_meta_L5.txt',
                  'data/tax_16S_uparse/prop_tables/fiveuL1201_demuxed_tax_normalized_otus_merged_meta_L6.txt')

uparse.tax <- lapply(files.uparse, function(x) read.table(x, header=TRUE, sep='\t', row.names=1))

names(uparse.tax) <- c('phylum', 'class', 'order', 'family', 'genus')

## rename columns with master sample IDs
uparse.tax <- lapply(seq_along(uparse.tax), function(x){
  colnames(uparse.tax[[x]]) <- sample.key$master_sampleid[match(colnames(uparse.tax[[x]]),
                                                                sample.key$uparse_sampleid)]
  uparse.tax[[x]]})

names(uparse.tax) <- c('phylum', 'class', 'order', 'family', 'genus')

## reduce dataset to exclude '5 µL' samples
uparse.tax.red <- lapply(seq_along(uparse.tax), function(x){
  pick <- colnames(uparse.tax[[x]])[colnames(uparse.tax[[x]]) %in% sample.key.red$master_sampleid]
  uparse.tax[[x]][,pick]})

names(uparse.tax.red) <- c('phylum', 'class', 'order', 'family', 'genus')

## determine which sample IDs overlap between the two datasets (we can do this at one level)
overlap.16S.id <- colnames(uparse.tax.red$phylum)[colnames(uparse.tax.red$phylum) %in% colnames(qiime.tax.red$phylum)]

## determine which taxa overlap between the two datasets (do this for each level)
overlap.16S.tax <- lapply(seq_along(uparse.tax.red), function(x){
  rownames(uparse.tax.red[[x]])[rownames(uparse.tax.red[[x]]) %in% rownames(qiime.tax.red[[x]])]})

## reduce datasets for comparison (+ transpose so samples are in rows)
qiime.compare <- lapply(seq_along(qiime.tax.red), function(x){
  dat <- as.data.frame(t(qiime.tax.red[[x]][overlap.16S.tax[[x]], overlap.16S.id]))
  # dat$Other <- 1-rowSums(dat) ## uncomment this line to group excluded taxa as 'Other'
  dat})

uparse.compare <- lapply(seq_along(uparse.tax.red), function(x){
  dat <- as.data.frame(t(uparse.tax.red[[x]][overlap.16S.tax[[x]], overlap.16S.id]))
  # dat$Other <- 1-rowSums(dat) ## uncomment this line to group excluded taxa as 'Other'
  dat})

## compute NMDS and perform Procrustes rotation to compare datasets
# NMDS
qiime.compare.mds <- lapply(qiime.compare, metaMDS)
uparse.compare.mds <- lapply(uparse.compare, metaMDS)

## plot procrustes rotations (UPARSE rotated relative to QIIME)
dir.create('results/procrustes_out')
pdf('results/procrustes_out/nmds_uparse_vs_qiime.pdf', width = 8, height = 6, pointsize = 8)
par(mfrow=c(2,3))
lapply(seq_along(qiime.compare.mds), function(x){
  plot(procrustes(qiime.compare.mds[[x]], uparse.compare.mds[[x]]), 
       xlab = 'NMDS1', ylab = 'NMDS2',
       ar.col=c('purple','blue','green','orange','red')[x],
       main = paste('Procrustes rotation \n(', tax.key$level[x], ' level)', sep=''))})
dev.off()

## plot procrustes rotations (UPARSE rotated relative to QIIME), all plots on same axis scale
pdf('results/procrustes_out/nmds_uparse_vs_qiime_monoscale.pdf', width = 8, height = 6, pointsize = 8)
par(mfrow=c(2,3))
lapply(seq_along(qiime.compare.mds), function(x){
  plot(procrustes(qiime.compare.mds[[x]], uparse.compare.mds[[x]]), 
       xlab = 'NMDS1', ylab = 'NMDS2',
       ar.col=c('purple','blue','green','orange','red')[x],
       xlim = c(-3,1), ylim = c(-2,2),
       main = paste('Procrustes rotation \n(', tax.key$level[x], ' level)', sep=''))})
dev.off()

## save R object
save.image('results/procrustes_out/ger_procrustes.RData')
