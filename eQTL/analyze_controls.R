library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggpubr)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')

gc()
output_fp = 'outputs/TEST_BRCA'
output_dirs = list.files(output_fp)
coll = data.frame()
for (i in output_dirs){
  if (i == 'readme.txt'){next}
  if (i == 'tcga_luad_raw_seq'){next}
  trans_eqtls = read.table(file.path(output_fp, i, 'trans_eqtl_fdr005.txt'), sep='\t', header = T)
  trans_eqtls['experiment'] = i
  coll = rbind(coll, trans_eqtls)
}

meta_mut_dir = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_muts = load_meta_mut(file.path(meta_mut_dir, list.files(meta_mut_dir)))

library(xlsx)
ctrs_raw = read.xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx', sheetIndex = 1)
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]

for (epr in unique(coll$experiment)){
  sub_epr = subset(coll, experiment == epr)
  for (snp in unique(sub_epr$protein_change)){
    nm = paste(epr, snp, sep='|')
    ctrs[[paste(nm, 'FDR', sep='|')]] = NA
    ctrs[[paste(nm, 'beta', sep='|')]] = NA
    sub_c = subset(sub_epr, protein_change==snp)
    rownames(sub_c) = sub_c$gene
    hit = intersect(ctrs$Gene, sub_c$gene)
    for (h in hit){
      ctrs[which(ctrs$Gene==h), paste(nm, 'FDR', sep='|')] = sub_c[h, 'FDR']
      ctrs[which(ctrs$Gene==h), paste(nm, 'beta', sep='|')] = sub_c[h, 'beta']
    }
  }
}

ctrs_long = gather(ctrs, key='group', value='value', 3:ncol(ctrs))
ctr_long = cbind(ctrs_long, t(as.data.frame(strsplit(ctrs_long$group, split = '\\|'))))
rownames(ctr_long) = NULL
colnames(ctr_long)[(ncol(ctr_long)-2):ncol(ctr_long)] = c('experiment', 'mutation', 'type')
ctr_long = ctr_long[,-which(colnames(ctr_long)=='group')]
ctr_long = spread(ctr_long, key=type, value=value)

gene_order = ctrs$Gene[order(ctrs$Gene_annotation, ctrs$Gene)]
ctr_long$Gene = factor(ctr_long$Gene, gene_order)
ctr_long['wrong'] = NA
ctr_long[which(ctr_long$beta<0), 'wrong'] = T
ctr_long[which(ctr_long$beta>-0), 'wrong'] = F

g = ggplot(ctr_long) +
  geom_point(aes(x=Gene, y=mutation, size=abs(beta), 
                 color=-log2(FDR), shape=wrong)) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~experiment, ncol = 2) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  mytme +
  labs(x='', y='Mutation (group)') +
  theme(strip.background = element_rect(fill=NA), 
        strip.text = element_text(face='bold', color='black', size=12),
        axis.text.x = element_text(hjust = 1, size=12),
        axis.text.y = element_text(hjust = 1, size=12),
        axis.title.y = element_text(size=16),
        panel.border = element_rect(color='black', size=1, fill=NA))
g
ggsave('/Users/jefft/Desktop/p53_project/Plots/eQTL/test2.pdf',
       plot=g, width=14,height=8,units='in',device='pdf',dpi=300)

