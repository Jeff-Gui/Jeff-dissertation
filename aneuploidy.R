library(tidyverse)
library(ggpubr)
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
setwd('/Users/jefft/Desktop/p53_project/CCLE')
dt = read.csv('RawData/aneuploidy_scores.csv')
mut = read.csv('ProcessedData/cl_meta_p53mut_merged.csv')

pool = intersect(unique(mut$depmap_id), dt$DepMap_ID)
wt_cl = intersect(pool, unique(mut$depmap_id[which(is.na(mut$Variant.Type))]))
mut_cl = pool[which(!pool %in% wt_cl)]
dt['p53_mut'] = '-'
dt$p53_mut[which(dt$DepMap_ID %in% mut_cl)] = '+'
lin = mut$lineage_1
names(lin) = mut$depmap_id
dt['lineage'] = lin[dt$DepMap_ID]

dt_plt = subset(dt, dt$lineage=='Blood')
dt_plt = na.omit(dt)
ggplot(dt_plt, aes(x=p53_mut, y=Aneuploidy.score)) + mytme +
  geom_boxplot(width=0.5, outlier.colour = NA) +
  geom_point(position = position_jitter(width=0.1), size=0.5) +
  labs(x='p53 Mutation', y='Aneuploidy Score') +
  facet_wrap(~lineage, nrow=4) +
  stat_compare_means(label.y.npc = 0.9, label.x.npc = 0.5, 
                     label = 'p.signif')
ggsave('plots/aneu_mut.pdf', dpi=300, height = 10, width = 11)

