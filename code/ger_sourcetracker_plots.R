## Summarize Gerlinger SourceTracker results
## Roxana Hickey <roxana.hickey@gmail.com
## 2015-11-09

setwd('~/Documents/projects/gerlinger/')
library(ggplot2)
library(reshape)
library(RColorBrewer)

## input sourcetracker results (R data)
load('results/sourcetracker_out/ger_wood_barberan_genus_subsample50/results.RData')
res <- results
rm(results)

## input sourcetracker output map files
map.out <- read.table('results/sourcetracker_out/ger_wood_barberan_genus_subsample50/map.txt',
                      header = T, sep = '\t', row.names = 1, comment.char = '')

## note control samples
control.id <- rownames(map.out)[grep('oneuLcontroL', rownames(map.out))]

## reduce maps to source/sink proportions
prop.pick <- c('Env', 'TITLE', 'Proportion_Human.Feces', 'Proportion_Human.Mouth', 'Proportion_Human.Skin', 
               'Proportion_Human.Urine', 'Proportion_Human.Vagina', 'Proportion_Unknown')

res.prop <- map.out[map.out$SourceSink == 'sink', prop.pick]
res.prop$SampleID <- rownames(res.prop)

## melt dataframe for barplots
res.prop.lg <- melt(res.prop, id.vars = c('SampleID', 'Env', 'TITLE'))

## narrow selection to Gerlinger only
res.prop.lg.ger <- subset(res.prop.lg, TITLE == 'Gerlinger')
sample.pick.ger <- rownames(map.out)[grep('oneuL', rownames(map.out))]

## note control samples
control.id <- rownames(map.out)[grep('oneuLcontroL', rownames(map.out))]

## set sample order for plots
ord.id <- c(rownames(map.out)[order(map.out$Env)][rownames(map.out)[order(map.out$Env)] %in% sample.pick.ger[!(sample.pick.ger %in% control.id)]],
            control.id)
ord.env <- c(as.character(map.out$Env[order(map.out$Env)][rownames(map.out)[order(map.out$Env)] %in% sample.pick.ger[!(sample.pick.ger %in% control.id)]]),
             rep('control', 2))

## bar plot of Gerlinger-only sourcetracker results
gg.prop.bar <- ggplot(res.prop.lg.ger, aes(x = SampleID, y = value, fill = variable))
gg.prop.bar + geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(brewer.pal(11, 'Spectral')[c(11,10,9,8,7)], 
                               'grey90'),
                    labels = c('Human Feces', 'Human Mouth', 'Human Skin',
                               'Human Urine', 'Human Vagina', 'Unknown'),
                    name = 'Putative Source') +
  scale_x_discrete(limits = rev(ord.id), 
                   labels = rev(ord.env)) +
  xlab('Room Type') +
  ylab('Proportion') +
  coord_flip() +
  theme_bw()
ggsave('figures/ger_sourcetracker_genus.png', width = 8, height = 10)

## summarize all human together
res.prop.2 <- res.prop
res.prop.2$Proportion_Human <- rowSums(res.prop.2[ , c('Proportion_Human.Feces',
                                                       'Proportion_Human.Mouth',
                                                       'Proportion_Human.Skin',
                                                       'Proportion_Human.Urine',
                                                       'Proportion_Human.Vagina')])
res.prop.2 <- res.prop.2[ , c('Env', 'TITLE', 'SampleID', 
                              'Proportion_Human', 'Proportion_Unknown')]

## density plot of 3 studies
gg.prop.2.hist <- ggplot(res.prop.2, aes(x = Proportion_Human, fill = TITLE))
gg.prop.2.hist + geom_density(alpha = 0.6, adjust = 0.75) +
  # geom_line(stat = 'density', adjust = 0.75) +
  scale_fill_manual(values = c('darkgoldenrod4', 'yellow2', 'cadetblue2'),
                    limits = c('Gerlinger',                                                               
                               'Athletic equipment microbiota are shaped by interactions with human skin',
                               'The ecology of microscopic life in household dust'), 
                    labels = c('Gerlinger',
                               'Athletic Facilities',
                               'Homes'),
                    name = 'Study') +
  xlim(c(0,1)) +
  xlab('Proportion of putative human-sourced bacteria') + 
  ylab('Frequency') +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme_bw() +
  theme(legend.key = element_rect(colour = 'black'), legend.background = element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))
ggsave('figures/ger_wood_barberan_sourcetracker_genus_density.png', width = 6, height = 6)

## violin plot of 3 studies
gg.prop.2.vio <- ggplot(res.prop.2, aes(x = TITLE, y = Proportion_Human, fill = TITLE))
gg.prop.2.vio + geom_violin(adjust = 0.75) +
  geom_boxplot(width = 0.1, fill = 'black', outlier.color = NA) +
  stat_summary(fun.y = median, geom = 'point', fill = 'white', shape = 21, size = 2.5) +
  scale_fill_manual(values = c('darkgoldenrod4', 'yellow2', 'cadetblue2'),
                    limits = c('Gerlinger',                                                               
                               'Athletic equipment microbiota are shaped by interactions with human skin',
                               'The ecology of microscopic life in household dust'), 
                    labels = c('Gerlinger',
                               'Athletic Facilities',
                               'Homes'),
                    name = 'Study') +
  scale_x_discrete(limits = c('Gerlinger',                                                               
                              'Athletic equipment microbiota are shaped by interactions with human skin',
                              'The ecology of microscopic life in household dust'), 
                   labels = c('Gerlinger',
                              'Athletic Facilities',
                              'Homes')) +
  xlab('Proportion of putative human-sourced bacteria') + 
  ylab('Frequency') +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme_bw() +
  theme(legend.key = element_rect(colour = 'black'), legend.background = element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))
ggsave('figures/ger_wood_barberan_sourcetracker_genus_violin.png', width = 6, height = 6)


## box plot of 3 studies
# gg.prop.2.box <- ggplot(res.prop.2, aes(x = TITLE, y = Proportion_Human, fill = TITLE))
# gg.prop.2.box + geom_boxplot(show_guide = FALSE) + 
#   scale_fill_manual(values = c('darkgoldenrod4', 'yellow2', 'cadetblue2'),
#                     limits = c('Gerlinger',                                                               
#                                'Athletic equipment microbiota are shaped by interactions with human skin',
#                                'The ecology of microscopic life in household dust'), 
#                     labels = c('Gerlinger',
#                                'Athletic Facilities',
#                                'Homes'),
#                     name = 'Study') +
#   scale_x_discrete(limits = c('Gerlinger',                                                               
#                               'Athletic equipment microbiota are shaped by interactions with human skin',
#                               'The ecology of microscopic life in household dust'), 
#                    labels = c('Gerlinger',
#                               'Athletic Facilities',
#                               'Homes')) +
#   xlab('Study') +
#   ylab('Proportion of putative human-sourced bacteria') +
#   theme_bw()
# ggsave('figures/ger_wood_barberan_sourcetracker_genus_boxplot.png', width = 6, height = 8)