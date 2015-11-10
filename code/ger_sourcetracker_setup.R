## Set up Gerlinger OTU data
## Roxana Hickey <roxana.hickey@gmail.com
## 2015-10-14

# setwd('~/Documents/projects/gerlinger/')

## read in mapping file (Gerlinger + HMP)
map.ger <- read.table('data/sourcetracker/merge_mapping_files/151020_ger_sourcetracker_map_noindoor.txt',
                      header = T, sep = '\t', row.names = 1, comment.char = '')

## read in OTU data for GENUS (L6) level**
## ** Note: this is filtered to OTUs that occur only in 2 or more samples
otu.ger.genus <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L6_filtered_2samples.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

otu.hmp.genus <- read.table('data/sourcetracker/sourcetracker_tutorial_files/hmp_otu_taxonomy_summaries/hmp_otu_table_L6_nospaces.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

length(rownames(otu.ger.genus)[rownames(otu.ger.genus) %in% rownames(otu.hmp.genus)]) # 254
length(rownames(otu.hmp.genus)[rownames(otu.hmp.genus) %in% rownames(otu.ger.genus)]) # 254

## merge OTU data for Gerlinger and HMP
otu.merge.genus <- merge(otu.ger.genus, otu.hmp.genus, by = 'row.names', all = TRUE)
row.names(otu.merge.genus) <- otu.merge.genus$Row.names

## reduce OTU table to samples present in mapping table
otu.merge.genus <- otu.merge.genus[ ,(colnames(otu.merge.genus) %in% row.names(map.ger))]

## remove contaminant OTUs (chloroplasts, mitochondria, Halomonas**)
## ** based on analysis of Gerlinger data; see ger_contaminants.Rmd
plants <- grep('Chloroplast', row.names(otu.merge.genus), ignore.case = TRUE)
mito <- grep('mitochondria', row.names(otu.merge.genus), ignore.case = TRUE)
halo <- grep('Halomonadaceae', row.names(otu.merge.genus), ignore.case = TRUE)
cont <- c(halo, mito, plants)

otu.merge.genus <- otu.merge.genus[-cont, ]

## write merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change 
## the first cell of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.genus, 'data/sourcetracker/merge_hmp_ger_16S_uparse/151109_merged_ger_hmp_L6_filtered_2samples.txt', 
            sep = '\t', quote = FALSE)

####################################
## Combine merged Gerlinger/Wood/Barberan data with HMP source

## read in merged OTU table (L6, genus)
## ** Note: this includes the Gerlinger samples filtered to OTUs that occur only in 2+ samples,
## athletic gym data (Wood et al. 2015), and 1500 homes data (Barberan et al. 2015). These three
## studies sequenced the 16S rRNA V4 region and used UPARSE for OTU picking + Greengenes for classification.
otu.ger.wood.barb <- read.table('data/sourcetracker/merge_all_datasets/151109_16S_otu_merge_ger_wood_barberan.txt', 
                                header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')
colnames(otu.ger.wood.barb) <- gsub('X', '', colnames(otu.ger.wood.barb))

## read in merged mapping file
map.ger.wood.barb <- read.table('data/sourcetracker/merge_mapping_files/151109_ger_barberan_wood_hmp_merged_sourcetracker_map.txt',
                                header = T, sep = '\t', row.names = 1, comment.char = '')

length(rownames(otu.ger.wood.barb)[rownames(otu.ger.wood.barb) %in% rownames(otu.hmp.genus)]) # 412
length(rownames(otu.hmp.genus)[rownames(otu.hmp.genus) %in% rownames(otu.ger.wood.barb)]) # 412

## merge OTU data for Gerlinger and HMP
otu.merge.all <- merge(otu.ger.wood.barb, otu.hmp.genus, by = 'row.names', all.x = TRUE, all.y = TRUE)
row.names(otu.merge.all) <- otu.merge.all$Row.names

## reduce OTU table to samples present in mapping table
otu.merge.all <- otu.merge.all[ ,(colnames(otu.merge.all) %in% row.names(map.ger.wood.barb))]

## remove contaminant OTUs (chloroplasts, mitochondria, Halomonas**)
## ** based on analysis of Gerlinger data; see ger_contaminants.Rmd
plants <- grep('Chloroplast', row.names(otu.merge.all), ignore.case = TRUE)
mito <- grep('mitochondria', row.names(otu.merge.all), ignore.case = TRUE)
halo <- grep('Halomonadaceae', row.names(otu.merge.all), ignore.case = TRUE)
cont <- c(halo, mito, plants)

## exclude contaminants, samples with fewer than 400 sequences, and OTUs with 0 hits
otu.merge.all <- otu.merge.all[-cont, ]
otu.merge.all <- otu.merge.all[rowSums(otu.merge.all, na.rm = T) > 0, 
                               colSums(otu.merge.all, na.rm = T) >= 400]

## write merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change 
## the first cell of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.all, 'data/sourcetracker/merge_all_datasets/151109_merged_ger_filtered_2samples_wood_barberan_hmp_L6.txt', 
            sep = '\t', quote = FALSE)

####################################
## Generate random subsample of Wood and Barberan studies to run shorter SourceTracker analysis
table(map.ger.wood.barb$TITLE) # view all projects in map table

pick.ger <- which(map.ger.wood.barb$TITLE == 'Gerlinger')

pick.hmp <- which(map.ger.wood.barb$TITLE == 'Bacterial Community Variation in Human Body Habitats Across Space and Time')

pick.wood <- which(map.ger.wood.barb$TITLE == 'Athletic equipment microbiota are shaped by interactions with human skin')
samp.wood <- sample(pick.wood, 50)

pick.barberan <- which(map.ger.wood.barb$TITLE == 'The ecology of microscopic life in household dust')
samp.barberan <- sample(pick.barberan, 50)

map.ger.wood.barb.subsamp50 <- map.ger.wood.barb[c(pick.hmp, pick.ger, samp.wood, samp.barberan), ]

write.table(map.ger.wood.barb.subsamp50,
            'data/sourcetracker/merge_mapping_files/151109_ger_barberan_wood_hmp_merged_sourcetracker_map_subsample50.txt',
            sep = '\t', quote = FALSE)

####################################
## Gerlinger/HMP OTU data for FAMILY (L5) level
# otu.ger.family <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L5.txt',
# header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

# otu.hmp.family <- read.table('data/sourcetracker_16S_uparse/sourcetracker_tutorial_files/hmp_otu_taxonomy_summaries/hmp_otu_table_L5_nospaces.txt',
# header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

# length(rownames(otu.ger.family)[rownames(otu.ger.family) %in% rownames(otu.hmp.family)]) # 125
# length(rownames(otu.hmp.family)[rownames(otu.hmp.family) %in% rownames(otu.ger.family)]) # 125

# otu.merge.family <- merge(otu.ger.family, otu.hmp.family, by = 'row.names', all = TRUE)

## FILTERED otu data for FAMILY (L5) level
# otu.ger.family.fil <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L5_filtered_2samples.txt',
# header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

# length(rownames(otu.ger.family.fil)[rownames(otu.ger.family.fil) %in% rownames(otu.hmp.family)]) # 114
# length(rownames(otu.hmp.family)[rownames(otu.hmp.family) %in% rownames(otu.ger.family.fil)]) # 114

# otu.merge.family.fil <- merge(otu.ger.family.fil, otu.hmp.family, by = 'row.names', all = TRUE)