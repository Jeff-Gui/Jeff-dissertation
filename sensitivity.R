library(tidyverse)
library(patchwork)
setwd('/Users/jefft/Desktop/p53_project')

dt = read.csv('CCLE/RNAi_TP53.csv', header = T, row.names = 1)
drug = read.csv('CCLE/drug_19Q4.csv', header = T, row.names = 1)
hist(dt$TP53, breaks = 100)

ggplot(dt) +
  theme_classic() +
  geom_boxplot(aes(x=TP53))

# Cell and drugs of interest
cell_set = c('MDAMB231', 'SUM149PT', 'HCC1395', 'MDAMB468', 'BT549')
drug_set = c('bortezomib', 'carfilzomib')

drug_idx = sapply(drug_set, function(x){
  return(grep(x, colnames(drug)))
})
cell_idx = sapply(cell_set, function(x){
  return(rownames(dt)[which(dt$cell_line_display_name==x)])
})
drug_OI = drug[,drug_idx]
colnames(drug_OI) = drug_set

co_rownm = rownames(drug_OI)[which(rownames(drug_OI) %in% rownames(dt))]

drug_and_RNAi = cbind(dt[co_rownm,], drug_OI[co_rownm,])
drug_and_RNAi['is_COI'] = FALSE
drug_and_RNAi[cell_idx,'is_COI'] = TRUE
a = ggplot(drug_and_RNAi, aes(x=TP53, y=-carfilzomib)) +
  theme_classic() +
  geom_point(aes(color=is_COI, alpha=is_COI), size=2) +
  geom_text(aes(label=cell_line_display_name, size=is_COI), nudge_x = 0.1, nudge_y = 0.17) +
  scale_size_manual(values = c(0,3)) +
  scale_alpha_manual(values=c(0.2, 1)) +
  theme(text = element_text(face='bold'), legend.position = 'none') +
  labs(title='CCLE', x='TP53 dependency score (RNAi)', y='Carfilzomib sensitivity (-log2FC)')
  #geom_smooth(method = 'lm')

b = ggplot(drug_and_RNAi, aes(x=TP53, y=-bortezomib)) +
  theme_classic() +
  geom_point(aes(color=is_COI, alpha=is_COI), size=2) +
  geom_text(aes(label=cell_line_display_name, size=is_COI), nudge_x = 0.1, nudge_y = 0.17) +
  scale_size_manual(values = c(0,3)) +
  scale_alpha_manual(values=c(0.2, 1)) +
  theme(text = element_text(face='bold'), legend.position = 'none') +
  labs(title='CCLE', x='TP53 dependency score (RNAi)', y='Bortezomib sensitivity (-log2FC)')

a+b

