library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')

gc()
output_fp = 'outputs/TEST_BRCA'
output_fp = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult/outputs'
coll = load_eQTL_output(output_fp, mode='fdr-no', exclude = 'tcga_nine_pool')
if (length(which(is.na(coll))>0)){
  coll = coll[-which(is.na(coll))]
}
coll = Reduce(rbind, coll)
# config which mutation group(s) to plot
table(coll$protein_change)
name_map = read.csv('/Users/jefft/Desktop/p53_project/datasets/mut_name_map.csv', header=T)
mutation_to_plt = name_map$label
names(mutation_to_plt) = name_map$code

mutation_to_plt = mutation_to_plt[c('hot_spot', 'contact', 'core', 'p.R273H', 'p.R175H', 'breast_w2016')]

all_mut = names(table(coll$protein_change))
mut_exclude = all_mut[which(!all_mut %in% name_map$code)]
print(mut_exclude)  # this mutant will not be ploted

meta_mut_dir = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_muts = load_meta_mut(file.path(meta_mut_dir, list.files(meta_mut_dir)))

library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
ctrs_clean = ctrs

ctrs = ctrs_clean
for (epr in unique(coll$experiment)){
  sub_epr = subset(coll, experiment == epr)
  for (snp in unique(names(mutation_to_plt))){
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

# coll$ann = NA
# coll$ann[coll$gene %in% unique(ctrs$Gene)] = 'hit'
# coll = subset(coll, coll$FDR < 0.05 & coll$experiment=='tcga_brca_raw_seq')
# ggplot() +
#   geom_point(data = subset(coll, is.na(coll$ann)), 
#              aes(x=beta, y=-log10(FDR)),size=0.1) +
#   geom_point(data = subset(coll, coll$ann == 'hit'),
#              aes(x=beta, y=-log10(FDR), color='red'), size=2)

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

# ctr_long = subset(ctr_long, ctr_long$mutation %in% names(mutation_to_plt))
ctr_long$mutation = mutation_to_plt[ctr_long$mutation]
ctr_long$mutation = factor(ctr_long$mutation, levels = mutation_to_plt)

# ctr_long$beta[which(ctr_long$beta <= 1)] = NA
plt.list = list()
# normalize beta to 0~1 range for vis comfort
# ctr_long$beta = 3 *(ctr_long$beta - min(ctr_long$beta, na.rm = T)) / (max(ctr_long$beta, na.rm = T) - min(ctr_long$beta, na.rm = T))
# myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette = colorRampPalette(c("royalblue","purple", "coral"))
up_p = ceiling(-log10(min(ctr_long$FDR, na.rm = T)))
print(up_p)
sc = scale_colour_gradientn(colours = myPalette(100), 
            limits=c(0, up_p),breaks = seq(0,12,4), trans='log1p')
for (i in unique(ctr_long$experiment)){
  ctr_long_sub = subset(ctr_long, ctr_long$experiment==i)
  g = ggplot(ctr_long_sub) +
    geom_point(aes(x=Gene, y=mutation, size=abs(beta), 
                   color=-log10(FDR), shape=wrong)) +
    scale_size_continuous(limits = c(0,ceiling(max(ctr_long$beta, na.rm=T)))) +
    # facet_rep_wrap(~experiment, ncol = 2, repeat.tick.labels = TRUE) +
    # scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
    sc +
    mytme +
    geom_point(data = subset(ctr_long_sub, FDR < 0.05), 
               aes(x=Gene, y=mutation), shape = '*', size=4, color='palegreen') +
    # y = 'Mutation (group)'
    labs(x='', y='', title = i) +
    theme(strip.background = element_rect(fill=NA), 
          strip.text = element_text(face='bold', color='black', size=12),
          axis.text.x = element_text(hjust = 1, vjust=0.5, size=12, angle = 90, face='italic'),
          axis.text.y = element_text(hjust = 1, size=12),
          axis.title.y = element_text(size=16),
          panel.border = element_rect(color='black', size=1, fill=NA),
          panel.grid.major = element_line(color='grey',linetype='dotted'),
          legend.position = 'none'
          )
  plt.list[[i]] = g
}

marrangeGrob(grobs=plt.list,ncol=3,nrow=3) %>% 
  ggsave('/Users/jefft/Desktop/p53_project/Plots/eQTL/controls/TCGA-pan_VS-mutneg_ult-noFDR_test.pdf',
       plot=., width=11.69*2.4,height=8.27*1.5,units='in',device='pdf',dpi=300)
