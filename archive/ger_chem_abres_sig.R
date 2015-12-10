library(ggplot2)
library(reshape)

setwd('~/Dropbox (BioBE)/BioBE Team Folder/Gerlinger/')

## read in data
dat <- read.table('chems_abr_forR.txt', header = T, sep = '\t')
dat$X.ABR. <- NULL

## make long format
dat.lg.1 <- melt(dat, id.vars = c('samples', 'TCSavg', 'TCCavg', 'MePBavg', 
                                  'EtPBavg', 'PrPBavg', 'BuPBavg'))
colnames(dat.lg.1) <- c('samples', 'Triclosan', 'Triclocarban', 'Methylparaben',
                        'Ethylparaben', 'Propylparaben', 'Butylparaben',
                        'gene', 'rpkm')

dat.lg.2 <- melt(dat.lg.1, id.vars = c('samples', 'gene', 'rpkm'))
colnames(dat.lg.2) <- c(colnames(dat.lg.1[ , c(1,8:9)]), 'chem', 'conc')

## plot
theme_set(theme_bw())

# makeTransparent<-function(someColor, alpha=100)
# {
#   newColor<-col2rgb(someColor)
#   apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
#                                               blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
# }
# 
# makeTransparent('darkmagenta', alpha = 60)
# [1] "#8B008B3C"

gg.dat <- ggplot(dat.lg.2, aes(x = conc, y = rpkm, color = gene))
gg.dat + geom_point(pch = 17, size = 3) +
  facet_wrap(~ chem, nrow = 3, scales = 'free') +
  stat_smooth(method = 'lm', se = FALSE, show_guide = FALSE) +
  scale_color_manual(values = c('darkmagenta', '#8B008B3C'),
                     name = 'Gene Family') +
#   scale_alpha_discrete(range = c(0.5, 1),
#                        name = 'Gene Family') +
  ylim(c(0,50)) +
  xlab('Mean chemical concentration (ng/g)') +
  ylab('Relative abundance of reads mapped to gene (RPKM)') + 
  theme(legend.position = 'bottom', 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave('~/Desktop/test.png', width = 8, height = 10, units = 'in')
  