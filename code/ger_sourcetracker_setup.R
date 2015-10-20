## Set up Gerlinger OTU data
## Roxana Hickey <roxana.hickey@gmail.com
## 2015-10-14

# setwd('~/Documents/gerlinger/')

## read in master sample ID key
sample.key <- read.table('data/sampleid_master_key.txt', header = T, sep = '\t')

## exclude the "5 µL" samples (poor sample quality)
sample.key.red <- sample.key[grep('05ul', sample.key$master_sampleid, invert = TRUE),]

## sourcetracker mapping file
map <- read.table('data/sourcetracker_16S_uparse/ger_sourcetracker_map_151013.txt',
                  header = T, sep = '\t', comment.char = '')

## otu data for GENUS (L6) level
otu.ger.genus <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L6.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

otu.hmp.genus <- read.table('data/sourcetracker_16S_uparse/sourcetracker_tutorial_files/hmp_otu_taxonomy_summaries/hmp_otu_table_L6_nospaces.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

length(rownames(otu.ger.genus)[rownames(otu.ger.genus) %in% rownames(otu.hmp.genus)]) # 285
length(rownames(otu.hmp.genus)[rownames(otu.hmp.genus) %in% rownames(otu.ger.genus)]) # 285

otu.merge.genus <- merge(otu.ger.genus, otu.hmp.genus, by = 'row.names', all = TRUE)

## write merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change the first cell 
## of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.genus, 'data/sourcetracker_16S_uparse/ger_sourcetracker_merged/merged_ger_hmp_L6_manual_nospaces.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

## otu data for FAMILY (L5) level
otu.ger.family <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L5.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

otu.hmp.family <- read.table('data/sourcetracker_16S_uparse/sourcetracker_tutorial_files/hmp_otu_taxonomy_summaries/hmp_otu_table_L5_nospaces.txt',
                            header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

length(rownames(otu.ger.family)[rownames(otu.ger.family) %in% rownames(otu.hmp.family)]) # 125
length(rownames(otu.hmp.family)[rownames(otu.hmp.family) %in% rownames(otu.ger.family)]) # 125

otu.merge.family <- merge(otu.ger.family, otu.hmp.family, by = 'row.names', all = TRUE)

## write merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change the first cell 
## of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.family, 'data/sourcetracker_16S_uparse/ger_sourcetracker_merged/merged_ger_hmp_L5_manual_nospaces.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

####################################
## now do this for the OTU tables filtered to OTUs that occur only in 2 or more samples
## FILTERED otu data for GENUS (L6) level
otu.ger.genus.fil <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L6_filtered_2samples.txt',
                                header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

length(rownames(otu.ger.genus.fil)[rownames(otu.ger.genus.fil) %in% rownames(otu.hmp.genus)]) # 254
length(rownames(otu.hmp.genus)[rownames(otu.hmp.genus) %in% rownames(otu.ger.genus.fil)]) # 254

otu.merge.genus.fil <- merge(otu.ger.genus.fil, otu.hmp.genus, by = 'row.names', all = TRUE)

## write FILTERED merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change the first cell 
## of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.genus.fil, 'data/sourcetracker_16S_uparse/ger_sourcetracker_merged/merged_ger_hmp_L6_manual_nospaces_filtered_2samples.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

## FILTERED otu data for FAMILY (L5) level
otu.ger.family.fil <- read.table('data/tax_16S_uparse/count_tables/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L5_filtered_2samples.txt',
                                 header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

length(rownames(otu.ger.family.fil)[rownames(otu.ger.family.fil) %in% rownames(otu.hmp.family)]) # 114
length(rownames(otu.hmp.family)[rownames(otu.hmp.family) %in% rownames(otu.ger.family.fil)]) # 114

otu.merge.family.fil <- merge(otu.ger.family.fil, otu.hmp.family, by = 'row.names', all = TRUE)

## write FILTERED merged OTU table
## *** NOTE *** you need to edit the output file to add a comment on the first line, and change the first cell 
## of the table to '#OTU ID'. I also replaced NAs with zeros in the output table.
write.table(otu.merge.family.fil, 'data/sourcetracker_16S_uparse/ger_sourcetracker_merged/merged_ger_hmp_L5_manual_nospaces_filtered_2samples.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)
