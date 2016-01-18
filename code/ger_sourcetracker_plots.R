## Summarize SourceTracker results
## Roxana Hickey <roxana.hickey@gmail.com
## Last updated 2015-12-09

setwd('~/Documents/projects/dust_2015/')
library(ggplot2)
library(reshape)
library(RColorBrewer)

## input sourcetracker results (R data)
# load('results/sourcetracker_out/ger_wood_barberan_genus_subsample50/results.RData')
load('results/sourcetracker_out/ger_wood_barberan_genus/results.RData')
res <- results
rm(results)

## input OTU/map table (to use map table)
load('results/otu_setup/ger_rm_contaminants.RData')

## input sourcetracker output map files
# map.out <- read.table('results/sourcetracker_out/ger_wood_barberan_genus_subsample50/map.txt',
map.out <- read.table('results/sourcetracker_out/ger_wood_barberan_genus/map.txt',
                      header = T, sep = '\t', row.names = 1, comment.char = '')
ger.map$Env <- gsub('.', ' ', ger.map$SpaceType, fixed = TRUE)
map.out$SpaceTypeBioBE <- ger.map$SpaceTypeBioBE[match(map.out$Env, ger.map$Env)]

## exclude controls and sample 1201 (low read count)
sample.ignore <- c('oneuL1201', 'oneuLcontroL1', 'oneuLcontroL2')
map.out.orig <- map.out
map.out <- map.out.orig[!(rownames(map.out.orig) %in% sample.ignore),]

## note control samples
# control.id <- rownames(map.out)[grep('oneuLcontroL', rownames(map.out))]

## reduce maps to source/sink proportions
prop.pick <- c('Env', 'SpaceTypeBioBE', 'TITLE', 'Proportion_Human.Feces', 'Proportion_Human.Mouth', 'Proportion_Human.Skin', 
               'Proportion_Human.Urine', 'Proportion_Human.Vagina', 'Proportion_Unknown')

res.prop <- map.out[which(map.out$SourceSink == 'sink'), prop.pick]
res.prop$SampleID <- rownames(res.prop)

## melt dataframe for barplots
res.prop.lg <- melt(res.prop, id.vars = c('SampleID', 'Env', 'SpaceTypeBioBE', 'TITLE'))

## narrow selection to present study only
res.prop.lg.ger <- subset(res.prop.lg, TITLE == 'University mixed-use facility')
sample.pick.ger <- rownames(map.out)[grep('oneuL', rownames(map.out))]

## note control samples
# control.id <- rownames(map.out)[grep('oneuLcontroL', rownames(map.out))]

## set sample order for plots
# ord.id <- c(rownames(map.out)[order(map.out$SpaceTypeBioBE)][rownames(map.out)[order(map.out$SpaceTypeBioBE)] %in% sample.pick.ger[!(sample.pick.ger %in% control.id)]],
#             control.id)
# ord.env <- c(as.character(map.out$SpaceTypeBioBE[order(map.out$SpaceTypeBioBE)][rownames(map.out)[order(map.out$SpaceTypeBioBE)] %in% sample.pick.ger[!(sample.pick.ger %in% control.id)]]),
#              rep('control', 2))
ord.id <- rownames(map.out)[order(map.out$SpaceTypeBioBE)][rownames(map.out)[order(map.out$SpaceTypeBioBE)] %in% sample.pick.ger]
ord.env <- as.character(map.out$SpaceTypeBioBE[order(map.out$SpaceTypeBioBE)][rownames(map.out)[order(map.out$SpaceTypeBioBE)] %in% sample.pick.ger])

## bar plot of sourcetracker results (present study only)
gg.prop.bar <- ggplot(res.prop.lg.ger, aes(x = SampleID, y = value, fill = variable))
gg.prop.bar + geom_bar(stat = 'identity') +
  # scale_fill_manual(values = c(brewer.pal(11, 'Spectral')[c(11,10,9,8,7)], 
                               # 'gainsboro'),
  scale_fill_manual(values = c(brewer.pal(9, 'BuPu')[c(9,8,7,6,5)], 
                               'gainsboro'),
                    labels = c('Human Feces', 'Human Mouth', 'Human Skin', 
                               'Human Vagina', 'Human Urine', 'Environmental/Other'),
                    name = 'Putative Source') +
  scale_x_discrete(limits = rev(ord.id), 
                   labels = rev(ord.env)) +
  xlab('Space Type') +
  ylab('Proportion') +
  coord_flip() +
  theme_bw()
ggsave('figures/ger_sourcetracker_genus_fulldata.png', width = 8, height = 10)
# ggsave('figures/ger_sourcetracker_genus.png', width = 8, height = 10)

## summarize all human together
res.prop.2 <- res.prop
res.prop.2$Proportion_Human <- rowSums(res.prop.2[ , c('Proportion_Human.Feces',
                                                       'Proportion_Human.Mouth',
                                                       'Proportion_Human.Skin',
                                                       'Proportion_Human.Urine',
                                                       'Proportion_Human.Vagina')])
res.prop.2 <- res.prop.2[ , c('Env', 'TITLE', 'SampleID', 
                              'Proportion_Human', 'Proportion_Unknown')]

## violin plot of 3 studies
gg.prop.2.vio <- ggplot(res.prop.2, aes(x = TITLE, y = Proportion_Human, color = TITLE, fill = TITLE))
gg.prop.2.vio + geom_violin(adjust = 0.75, alpha = 0.3) +
  geom_boxplot(width = 0.1, aes(fill = TITLE), color = 'black', outlier.color = NA) +
  stat_summary(fun.y = median, geom = 'point', fill = 'white', shape = 21, size = 2.5) +
  scale_fill_manual(values = c('darkmagenta', '#FF3300', '#FF9500'),
                    limits = c('University mixed-use facility',                                                               
                               'Athletic equipment microbiota are shaped by interactions with human skin',
                               'The ecology of microscopic life in household dust')) +
  scale_color_manual(values = c('darkmagenta', '#FF3300', '#FF9500'),
                    limits = c('University mixed-use facility',                                                               
                               'Athletic equipment microbiota are shaped by interactions with human skin',
                               'The ecology of microscopic life in household dust')) +
  scale_x_discrete(limits = c('University mixed-use facility',                                                               
                              'Athletic equipment microbiota are shaped by interactions with human skin',
                              'The ecology of microscopic life in household dust'), 
                   labels = c('Present Study',
                              'Athletic Facilities',
                              'Homes')) +
  ylab('Proportion of putative human-sourced bacteria') + 
  theme_bw() +
  theme(legend.position = 'none', axis.title.x = element_blank())
ggsave('figures/ger_wood_barberan_sourcetracker_genus_violin.png', width = 8, height = 8)

## box plot of 3 studies
# gg.prop.2.box <- ggplot(res.prop.2, aes(x = TITLE, y = Proportion_Human, fill = TITLE))
# gg.prop.2.box + geom_boxplot(show_guide = FALSE) + 
#   scale_fill_manual(values = c('darkmagenta', '#FFD701', 'turquoise'),
#                     limits = c('University mixed-use facility',                                                               
#                                'Athletic equipment microbiota are shaped by interactions with human skin',
#                                'The ecology of microscopic life in household dust'), 
#                     labels = c('University mixed-use facility',
#                                'Athletic Facilities',
#                                'Homes'),
#                     name = 'Study') +
#   scale_x_discrete(limits = c('University mixed-use facility',                                                               
#                               'Athletic equipment microbiota are shaped by interactions with human skin',
#                               'The ecology of microscopic life in household dust'), 
#                    labels = c('Present Study',
#                               'Athletic Facilities',
#                               'Homes')) +
#   xlab('Study') +
#   ylab('Proportion of putative human-sourced bacteria') +
#   theme_bw()
# ggsave('figures/ger_wood_barberan_sourcetracker_genus_boxplot.png', width = 6, height = 8)

## density plot of 3 studies
# gg.prop.2.hist <- ggplot(res.prop.2, aes(x = Proportion_Human, fill = TITLE))
# gg.prop.2.hist + geom_density(alpha = 0.6, adjust = 0.75) +
#   # geom_line(stat = 'density', adjust = 0.75) +
#   scale_fill_manual(values = c('darkmagenta', '#FFD701', 'turquoise'),
#                     limits = c('University mixed-use facility',                                                               
#                                'Athletic equipment microbiota are shaped by interactions with human skin',
#                                'The ecology of microscopic life in household dust'), 
#                     labels = c('Present Study',
#                                'Athletic Facilities',
#                                'Homes'),
#                     name = 'Study') +
#   xlim(c(0,1)) +
#   xlab('Proportion of putative human-sourced bacteria') + 
#   ylab('Frequency') +
#   guides(fill = guide_legend(override.aes = list(colour = NULL))) +
#   theme_bw() +
#   theme(legend.key = element_rect(colour = 'black'), legend.background = element_blank(),
#         legend.position = c(1,1), legend.justification = c(1,1))
# ggsave('figures/ger_wood_barberan_sourcetracker_genus_density_fulldata.png', width = 6, height = 6)
# ggsave('figures/ger_wood_barberan_sourcetracker_genus_density.png', width = 6, height = 6)