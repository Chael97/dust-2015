## Gerlinger OTU data
## Roxana Hickey <roxana.hickey@gmail.com
## 2015-10-14

# setwd('~/Documents/gerlinger/')

## read in master sample ID key
sample.key <- read.table('data/sampleid_master_key.txt', header = T, sep = '\t')

## exclude the "5 µL" samples (poor sample quality)
sample.key.red <- sample.key[grep('05ul', sample.key$master_sampleid, invert = TRUE),]

## sourcetracker mapping file
map <- read.table('data/sourcetracker_16S_uparse/ger_microbialsource_mappingfile_simplified_151013.txt',
                       header = T, sep = '\t', comment.char = '')

otu.ger.genus <- read.table('data/16S_uparse_otu_taxonomy/gerlinger_uparse_merged_taxonomy_summaries/gerlinger_uparse_merged_L6.txt',
                           header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

otu.hmp.genus <- read.table('data/sourcetracker_16S_uparse/sourcetracker_tutorial_files/hmp_otu_taxonomy_summaries/hmp_otu_table_L6_nospaces.txt',
                           header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

map$X.SampleID[!(map.test$X.SampleID %in% colnames(otu.merge.test))] # show sample names missing from OTU table

length(rownames(otu.ger.test)[rownames(otu.ger.test) %in% rownames(otu.hmp.test)])
length(rownames(otu.hmp.test)[rownames(otu.hmp.test) %in% rownames(otu.ger.test)])

otu.merge.test <- merge(otu.ger.test, otu.hmp.test, by = 'row.names', all = TRUE)
write.table(otu.merge.test, 'data/sourcetracker_16S_uparse/merged_ger_hmp_L6_manual_nospaces.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

# Note: you need to edit the output file to add a comment on the first line, and change the first cell of the table to '#OTU ID'
# I also replaced NAs with zeros in the output table.
