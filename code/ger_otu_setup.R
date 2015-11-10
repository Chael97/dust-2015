## Set up OTU and map tables for downstream analysis
## 2015-10-27
## Roxana Hickey <roxana.hickey@gmail.com>

setwd('~/Documents/projects/gerlinger/')
library(vegan)

## read in the genus-level OTU table
otu.g <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L6_filtered_2samples.txt',
                    header = TRUE, sep = '\t', row.names = 1, comment.char = '', skip = 1)

## read in mapping file with metadata
ger.map.orig <- read.table('data/map_files/151028_map_edit_finalchems_crack_plus.txt',
                           header = TRUE, sep = '\t', row.names = 1, comment.char = '')
ger.map.orig$SampleID <- rownames(ger.map.orig)

## reduce mapping file to oneuL samples
ger.map <- ger.map.orig[grep('oneuL', rownames(ger.map.orig)), ]

## define samples to ignore due to low read count or from another study
sample.ignore <- c('oneuL3031',
                   colnames(otu.g)[grep('oneuLP', colnames(otu.g))],
                   colnames(otu.g)[grep('oneuLC', colnames(otu.g))],
                   colnames(otu.g)[grep('oneuLFT', colnames(otu.g))])

## make reduced OTU table
otu.g <- otu.g[, !(colnames(otu.g) %in% sample.ignore)]

## make reduced mapping table
ger.map <- ger.map[!(rownames(ger.map) %in% sample.ignore),]

## look for obvious contaminants in controls to remove
cont.otu <- as.matrix(otu.g[, grep('controL', colnames(otu.g))])
cont.otu <- cont.otu[rowSums(cont.otu) > 0,] # exclude zero-sum OTUs
sort(rowSums(cont.otu), decreasing = TRUE)[1:10]
## the most abundant taxa by far are mitochondria (13,667 reads across both samples)
## and chloroplasts (13,439 reads), followed by Halomonadaceae (6,754 reads). 
## thereafter the read counts are below 2,000, or approximately 10% of each sample;
## control1 has 22,237 reads and control2 has 26,430 reads.
cont.exclude <- names(sort(rowSums(cont.otu), decreasing = TRUE))[1:3] # remove by order
cont.exclude <- rownames(cont.otu)[c(grep('chloroplast', rownames(cont.otu), ignore.case = TRUE),
                                     grep('mitochondria', rownames(cont.otu), ignore.case = TRUE),
                                     grep('halomonadaceae', rownames(cont.otu), ignore.case = TRUE))] # remove by name

## save original OTU table
otu.g.orig <- otu.g

## write new OTU table without contaminants
otu.g <- otu.g.orig[!(rownames(otu.g.orig) %in% cont.exclude),]

## make rarefied OTU count table
## first look at read counts across samples
summary(colSums(otu.g))

## rarefy to the minimum count of 6,452
otu.g.rare <- rrarefy(t(otu.g), 6400)

## transpose original OTU table to put samples in rows, OTUs in columns
otu.g <- t(otu.g)

## check for and remove zero-sum OTUs
zero.otus <- colnames(otu.g)[colSums(otu.g) == 0]
otu.g <- otu.g[, !(colnames(otu.g) %in% zero.otus)]

zero.otus.rare <- colnames(otu.g.rare)[colSums(otu.g.rare) == 0]
otu.g.rare <- otu.g.rare[, !(colnames(otu.g.rare) %in% zero.otus.rare)]

## replace rownames with shortened sample names (take out 'oneuL')
rownames(otu.g) <- ger.map$Description[match(rownames(otu.g), rownames(ger.map))]
rownames(otu.g.rare) <- ger.map$Description[match(rownames(otu.g.rare), rownames(ger.map))]
rownames(ger.map) <- ger.map$Description

## make proportion tables
otu.g.prop <- prop.table(otu.g, margin = 1)
otu.g.rare.prop <- prop.table(otu.g.rare, margin = 1)

## clean up unnecessary variables
rm(zero.otus, zero.otus.rare)

## save R data
# dir.create('results/otu_setup')
save.image('results/otu_setup/ger_otu_setup.RData')