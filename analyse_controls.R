setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult, TCGA-pan_VS-wt, TCGA-pan_VS-null
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'control')
data_out = file.path(dir_home, 'data_out')
library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')

gc()
# output_fp = 'outputs/TEST_BRCA'
output_fp = eqtl_out
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

mutation_to_plt = mutation_to_plt[c('hot_spot', 'contact', 'conformation', 'sandwich', 'p.R273H', 'p.R175H', 'breast_w2016')]

all_mut = names(table(coll$protein_change))
mut_exclude = all_mut[which(!all_mut %in% name_map$code)]
print(mut_exclude)  # this mutant will not be ploted

meta_mut_dir = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_muts = load_meta_mut(file.path(meta_mut_dir, list.files(meta_mut_dir)))

# load control genes
library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
# check if control genes are expressed
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
no_exp = c()
for (gn in unique(ctrs$Gene)){
  if (!gn %in% rownames(mtx)){
    no_exp = c(no_exp, gn)
  }
}
print(paste('Genes not detected:', paste(no_exp, collapse = ',')))
print(length(no_exp))
ctrs = ctrs[which(!ctrs$Gene %in% no_exp),]
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
ctr_long[which(ctr_long$beta>0), 'wrong'] = F

# ctr_long = subset(ctr_long, ctr_long$mutation %in% names(mutation_to_plt))
ctr_long$mutation = mutation_to_plt[ctr_long$mutation]
ctr_long$mutation = factor(ctr_long$mutation, levels = mutation_to_plt)

# ctr_long$beta[which(ctr_long$beta <= 1)] = NA
plt.list = list()
myPalette = colorRampPalette(c("royalblue","purple", "coral"))
up_p = ceiling(-log10(min(ctr_long$FDR, na.rm = T)))
print(up_p)
sc = scale_colour_gradientn(colours = myPalette(100), 
            limits=c(0, up_p),breaks = seq(0,12,4), trans='log1p')

ctr_long$cancer = toupper(sapply(ctr_long$experiment, function(x){return(strsplit(x, split='_')[[1]][2])}))

# see how many genes cen be recovered in each cancer, hotspot only
#### wt control one
ctr1 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control')
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==TRUE)),
                                        wrong_hit=length(which(wrong==FALSE)),
                                        not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
ftr_ctr1 = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_ctr1)
ctr1_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Negative control 1') +
  mytme
sm_tb_ctr1 = sm_tb[sm_tb$gp=='hit',]

#### wt control 2
ctr1 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control_2')
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==TRUE)),
                                                wrong_hit=length(which(wrong==FALSE)),
                                                not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
ftr_ctr2 = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_ctr2)
ctr2_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Negative control 2') +
  mytme
sm_tb_ctr2 = sm_tb[sm_tb$gp=='hit',]

#### positive controls
ctr1 = subset(ctr_long, !ctr_long$Gene_annotation %in% c('wt_control', 'wt_control_2'))
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==FALSE)),
                                                wrong_hit=length(which(wrong==TRUE)),
                                                not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
ftr_pos = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_pos)
pos_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Positive controls') +
  mytme
sm_tb_pos = sm_tb[sm_tb$gp=='hit',]

### Rank cancers
mtx = matrix(0, nrow=length(unique(sm_tb$cancer)), ncol=3)
rownames(mtx) = unique(sm_tb$cancer)
colnames(mtx) = c('ctr1', 'ctr2', 'pos')
mtx = as.data.frame(mtx)
for (i in 1:nrow(sm_tb_ctr1)){
  mtx[as.character(sm_tb_ctr1$cancer[i]), 'ctr1'] = sm_tb_ctr1$count[i]
}
for (i in 1:nrow(sm_tb_ctr2)){
  mtx[as.character(sm_tb_ctr2$cancer[i]), 'ctr2'] = sm_tb_ctr2$count[i]
}
for (i in 1:nrow(sm_tb_pos)){
  mtx[as.character(sm_tb_pos$cancer[i]), 'pos'] = sm_tb_pos$count[i]
}
mtx$sumHit = rowSums(mtx)
mtx = mtx[order(mtx$sumHit, decreasing = T),]
mtx$cancer = factor(rownames(mtx), levels=rownames(mtx))
rank_pfl = ggplot(mtx) +
  geom_bar(aes(x=cancer, y=sumHit), stat = 'identity', fill='black') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x='', y='Sum of hit counts') +
  mytme
marrangeGrob(grobs=list(ctr1_count, ctr2_count, pos_count, rank_pfl),ncol=2,nrow=2) %>% 
  ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_hitCount.pdf'),
         plot=., width=11.69*1.3,height=8.27,units='in',device='pdf',dpi=300)


# WT control or Not !!! Choice
ctr_long_plt = subset(ctr_long, !ctr_long$Gene_annotation %in% c('wt_control', 'wt_control_2'))
ctr_long_plt = subset(ctr_long, ctr_long$Gene_annotation=='wt_control')
ctr_long_plt = subset(ctr_long, ctr_long$Gene_annotation=='wt_control_2')

for (i in unique(ctr_long_plt$experiment)){
  ctr_long_sub = subset(ctr_long_plt, ctr_long_plt$experiment==i)
  g = ggplot(ctr_long_sub) +
    geom_point(aes(x=Gene, y=mutation, size=abs(beta), 
                   color=-log10(FDR), shape=wrong), stroke=2) +
    scale_shape_manual(name='', label=c('Negative', 'Positive'), limits = c(FALSE, TRUE), values = c(1,2)) +
    scale_size_continuous(limits = c(0,ceiling(max(abs(ctr_long_plt$beta), na.rm=T)))) +
    # facet_rep_wrap(~experiment, ncol = 2, repeat.tick.labels = TRUE) +
    # scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
    sc +
    mytme +
    geom_point(data = subset(ctr_long_sub, FDR < 0.05), 
               aes(x=Gene, y=mutation), shape = '*', size=6, color='black') +
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
  ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_posCtr.pdf'),
       plot=., width=11.69*2.4,height=8.27*1.5,units='in',device='pdf',dpi=300)


### Contact VS Conform in Esposito, 2022
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
load('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData')
tcga = dt
load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
ccle = dt
remove(dt)
gc()


contact_cell = intersect(ccle[[1]]@colData$NAME, 
                         c('MDA-MB-468','U373MG', 'U-251 MG', 'SF-295', 'HCC193', 'PC9'))
conform_cell = intersect(ccle[[1]]@colData$NAME, 
                         c('HCC1395', 'HCC1954', 'SK-MEL-2', 'SK-LMS-1'))
genes = c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4')

comp = data.frame('name'=c(contact_cell, conform_cell), 'group'=c(rep('contact', 3), rep('conform', 3)))
comp['ID'] = sapply(comp$name, function(x){ccle[[1]]@colData$PATIENT_ID[which(ccle[[1]]@colData$NAME==x)]})
comp = cbind(comp, t(assay(ccle[[1]][genes,comp$ID,'RNA'])))

comp = gather(comp, key='gene', value='exp', 4:7)
ggplot(comp, aes(x=group,y=exp)) +
  geom_point() +
  facet_wrap(~gene) +
  geom_text(aes(label=name)) +
  mytme

mutation_groups = list(c(273,248), c(175,245,249,282))
names(mutation_groups) = c('Contact', 'Conform')

genes = c('HMGCR', 'MVK', 'MVD', 'FDPS', 'SQLE', 'LSS', 'DHCR7', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4')
plt = get_genes_plt(genes=genes, ccle=ccle, tcga=tcga, mutation_groups, 
              primary_site = 'Breast', rnai = NULL, 
              comparison = list('rna'=list(c('Contact', 'Wildtype'),
                                           c('Conform', 'Contact')),
                                'rnai'=list(c('Contact', 'Wildtype'))),
              no_ccle = TRUE)
plt %>% marrangeGrob(ncol=2, nrow=3, top = '',
                     layout_matrix = matrix(1:6,byrow = T, ncol=2)) %>%
  ggsave(file.path(plot_out, 'TEAD_controls.pdf'),
         plot=., width=8,height=16,units='in',device='pdf',dpi=300)

