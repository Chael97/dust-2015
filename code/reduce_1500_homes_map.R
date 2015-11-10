## Reduce Barberan map table to only indoor samples for SourceTracker analysis
barb <- read.table('data/sourcetracker/barberan_1000homes/16S_otu_table_wTax/16S_otu_table_wTax_L6.txt',
                   header = T, sep = '\t', row.names = 1, skip = 1, comment.char = '')

barb.i <- colnames(barb)[grep('.I', colnames(barb))]
barb.i <- gsub('X', '', barb.i)

barb.map <- read.table('data/sourcetracker/barberan_1000homes/homes_mapping_file.txt',
                       header = T, sep = '\t', comment.char = '')
barb.map$SampleID <- paste0(barb.map$ID, '.I')

barb.map.indoor <- subset(barb.map, SampleID %in% barb.i)

write.table(barb.map.indoor, 'data/sourcetracker/barberan_1000homes/barberan_map_for_sourcetracker_indoor_only.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)