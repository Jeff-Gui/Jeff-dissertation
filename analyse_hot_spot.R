setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-wt, TCGA-pan_VS-mutneg_ult
# ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE/processed_dt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'hotspot')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../overlap_utils.R')
source('../revigo_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)
library(ggrepel)
library(ggsci)

# Load Gene ========
beta_cutoff = 0
coll = load_eQTL_output(eqtl_out, beta=beta_cutoff, exclude = c('tcga_nine_pool', 'metabric_raw_its'))
df_coll = data.frame()
for (i in names(coll)){
  nm = toupper(strsplit(i, split='_')[[1]][2])
  df = coll[[i]]
  if (class(df)!='logical'){
    df = subset(df, abs(df$beta) > beta_cutoff)
    df_coll = rbind(df_coll, df)
  }
}
df_coll[['cancer']] = sapply(df_coll$experiment, function(x){
  return(strsplit(x, split='_')[[1]][2])
})
df_coll$cancer = toupper(df_coll$cancer)

# Overlapping gene ========
## load control genes ====
library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
neg.1.ctrs = unique(ctrs[ctrs$Gene_annotation=='wt_control',]$Gene)
neg.2.ctrs = unique(ctrs[ctrs$Gene_annotation=='wt_control_2',]$Gene)
ctrs = ctrs[-grep('wt', ctrs$Gene_annotation),]
pos.ctrs = ctrs$Gene

library(ComplexUpset) # upset plot requires sample-group matrix
# https://blog.csdn.net/tuanzide5233/article/details/83109527
library(ggvenn)
#exps = unique(df_coll$experiment)
# exps = c('BLCA', 'STAD', 'BRCA', 'LGG', 'COAD')
exps = unique(df_coll$cancer)
hs = subset(df_coll, df_coll$protein_change=='hot_spot' & 
              df_coll$cancer %in% exps)
## Plot overview ====
hs_pos_smm = subset(hs, hs$beta > 0) %>% group_by(cancer) %>%
  summarise(count_pos = length(beta))
hs_neg_smm = subset(hs, hs$beta < 0) %>% group_by(cancer) %>%
  summarise(count_neg = length(beta))
hs_smm = full_join(hs_pos_smm, hs_neg_smm, by='cancer')
hs_smm[is.na(hs_smm)] = 0
hs_smm = rbind(hs_smm, c('LUAD',0,0))
hs_smm = rbind(hs_smm, c('OV',0,0))
hs_smm$count_pos = as.numeric(hs_smm$count_pos)
hs_smm$count_neg = as.numeric(hs_smm$count_neg)
hs_smm = hs_smm[order(hs_smm$count_pos,decreasing = T),]
hs_smm$cancer = factor(hs_smm$cancer, levels = hs_smm$cancer)
# normalize with mutant sample size
ssize = read.table(file.path(dir_home, 'sample_size.txt'), header = T)
ssize = ssize[ssize$Mutation=='hot_spot',]
hs_smm$mutant_size = 0
for (i in 1:nrow(hs_smm)){
  if (hs_smm$cancer[i] %in% ssize$Cancer){
    hs_smm$mutant_size[i] = ssize[ssize$Cancer==hs_smm$cancer[i],'Mutant']
  }
}
g = ggplot(hs_smm[!hs_smm$cancer %in% c('LUSC', 'LUAD', 'OV'),] %>% 
             gather(key='sign', value='Count', 2:3),
           aes(x=cancer,y=Count/mutant_size)) +
  geom_bar(aes(fill=sign),stat='identity', position = 'dodge') +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) +
  geom_text(aes(label=paste(Count, mutant_size, sep='/'), group=sign), 
            position = position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  labs(x='', y='Gene count / Mutant sample size') + scale_fill_d3(name='',labels = c('Negative', 'Positive')) +
  # geom_text(label=parse(text="'Not passing\nVAF filter'"), x='LUAD', y=2, size=2.5) +
  mytme + theme(legend.position = c(0.8,0.9))
ggsave(file.path(plot_out, 'panCan_count_overview.pdf'),
       plot=g, height=11.69*0.4,width=8.27*0.7,units='in',device='pdf',dpi=300)

# Identify common signature across cancers
coll_pos_hs = list()
for (i in exps){
  # !!! choose sign !!!
  sub_hs = subset(hs, hs$cancer==i & hs$beta > 0)
  if (nrow(sub_hs)>0 & !i %in% c('LUAD', 'LUSC')){
    sub_hs = sub_hs[order(abs(sub_hs$beta), decreasing = T),]
    # !!! choose which hit mode to use !!!
    # hits = sub_hs$gene
    hits = sub_hs$gene[1:min(1000,nrow(sub_hs))]
    coll_pos_hs[[i]] = hits
  }
}

# coll_pos_hs[['controls']] = pos.ctrs
upset_mtx_pan = gen_upSet_mtx(coll_pos_hs)
#upset_mtx_pan = upset_mtx_pan[,-which(colSums(upset_mtx_pan)==0)]
to_rm = which(rowSums(upset_mtx_pan)==0)
if (length(to_rm) > 0){
  upset_mtx_pan = upset_mtx_pan[-to_rm,]
}
print(dim(upset_mtx_pan))

### Test overlapping significance (dual) ====
set.seed(3.180110984)
colnames(upset_mtx_pan)
ovl_p = list()
for (i in 1:(ncol(upset_mtx_pan)-1)){
  for (j in (i+1):ncol(upset_mtx_pan)){
    tc = rep(NA, ncol(upset_mtx_pan))
    tc[i] = T
    tc[j] = T
    id = paste(colnames(upset_mtx_pan)[i], colnames(upset_mtx_pan)[j], sep=' X ')
    print(id)
    out = run_overlap_test(upset_mtx_pan, n_loop=1000,
                          test_condition = tc)
    ovl_p[[id]] = out
  }
}
save(ovl_p, file=file.path(data_out, 'overlap_test_Dual_panCan_hotspot_pos.RData'))

library(ggridges)
ridge_df = list()
for (i in 1:length(ovl_p)){
  ridge_df[[names(ovl_p)[i]]] = c(ovl_p[[i]]$null_dst, ovl_p[[i]]$observed)
}
ridge_df = as.data.frame(ridge_df)
colnames(ridge_df) = gsub('_', ' X ', colnames(ridge_df))
colnames(ridge_df) = gsub('\\.', ' ', colnames(ridge_df))
g = ggplot(ridge_df[-1,] %>% gather(key='gp',value='price',1:ncol(ridge_df)),
       aes(x=price, y=gp, fill=gp)) +
  geom_density_ridges() + labs(x='Count', y='') +
  geom_point(data = ridge_df[nrow(ridge_df),] %>% gather(key='gp', value='obs', 1:ncol(ridge_df)),
             aes(x=obs, y=gp), color='grey20', stroke=1) + mytme + 
  theme(legend.position = 'none', axis.text.y = element_text(size=10))
ggsave(file.path(plot_out, 'panCan_pos_oval_dualTest.pdf'),
       plot=g, height=11.69*0.5,width=8.27*0.5,units='in',device='pdf',dpi=300)

### Test overlapping significance (triple) ====
library(gtools)
set.seed(3.180110984)
cnm = colnames(upset_mtx_pan)
ovl_p = list()
cmb = combinations(ncol(upset_mtx_pan), 3)
for (i in 1:nrow(cmb)){
  tc = rep(NA, ncol(upset_mtx_pan))
  tc[c(cmb[i,1:3])] = T
  id = paste(cnm[cmb[i,1]], cnm[cmb[i,2]], cnm[cmb[i,3]], sep=' X ')
  print(id)
  out = run_overlap_test(upset_mtx_pan, n_loop=1000,
                         test_condition = tc)
  ovl_p[[id]] = out
}
save(ovl_p, file=file.path(data_out, 'overlap_test_Triple_panCan_hotspot_pos.RData'))

library(ggridges)
ridge_df = list()
for (i in 1:length(ovl_p)){
  ridge_df[[names(ovl_p)[i]]] = c(ovl_p[[i]]$null_dst, ovl_p[[i]]$observed)
}
ridge_df = as.data.frame(ridge_df)
colnames(ridge_df) = gsub('\\.', ' ', colnames(ridge_df))
g = ggplot(ridge_df[-1,] %>% gather(key='gp',value='price',1:ncol(ridge_df)),
           aes(x=price, y=gp, fill=gp)) +
  geom_density_ridges() + labs(x='Count', y='') +
  geom_point(data = ridge_df[nrow(ridge_df),] %>% gather(key='gp', value='obs', 1:ncol(ridge_df)),
             aes(x=obs, y=gp), color='grey20', stroke=1) + mytme + 
  theme(legend.position = 'none', axis.text.y = element_text(size=10))
ggsave(file.path(plot_out, 'panCan_pos_oval_tripleTest.pdf'),
       plot=g, height=11.69*0.8,width=8.27*0.8,units='in',device='pdf',dpi=300)


### Test overlapping significance (quard) ====
set.seed(3.180110984)
cnm = colnames(upset_mtx_pan)
ovl_p = list()
cmb = combinations(ncol(upset_mtx_pan), 4)
for (i in 1:nrow(cmb)){
  tc = rep(NA, ncol(upset_mtx_pan))
  tc[c(cmb[i,1:4])] = T
  id = paste(cnm[cmb[i,1]], cnm[cmb[i,2]], cnm[cmb[i,3]], cnm[cmb[i,4]], sep=' X ')
  print(id)
  out = run_overlap_test(upset_mtx_pan, n_loop=1000,
                         test_condition = tc)
  ovl_p[[id]] = out
}
save(ovl_p, file=file.path(data_out, 'overlap_test_Quard_panCan_hotspot_pos.RData'))

ridge_df = list()
for (i in 1:length(ovl_p)){
  ridge_df[[names(ovl_p)[i]]] = c(ovl_p[[i]]$null_dst, ovl_p[[i]]$observed)
}
ridge_df = as.data.frame(ridge_df)
colnames(ridge_df) = gsub('\\.', ' ', colnames(ridge_df))
g = ggplot(ridge_df[-1,] %>% gather(key='gp',value='price',1:ncol(ridge_df)),
           aes(x=price, y=gp, fill=gp)) +
  geom_density_ridges() + labs(x='Count', y='') +
  geom_point(data = ridge_df[nrow(ridge_df),] %>% gather(key='gp', value='obs', 1:ncol(ridge_df)),
             aes(x=obs, y=gp), color='grey20', stroke=1) + mytme + 
  theme(legend.position = 'none', axis.text.y = element_text(size=10))
ggsave(file.path(plot_out, 'panCan_pos_oval_quardTest.pdf'),
       plot=g, height=11.69*0.8,width=8.27*1.2,units='in',device='pdf',dpi=300)

### Plot upset ====
uquery = list()
count = 1
for (i in c('overlap_test_Quard_panCan_hotspot_pos.RData',
            'overlap_test_Triple_panCan_hotspot_pos.RData',
            'overlap_test_Dual_panCan_hotspot_pos.RData')){
  load(file.path(data_out, i))
  for (j in names(ovl_p)){
    if (ovl_p[[j]]$p.left >= 0.95){
      ist = strsplit(j, split=' X ')[[1]]
      uquery[[count]] = upset_query(intersect = ist, color='#FF7F0EFF', fill='#FF7F0EFF')
      count = count + 1
    }
  }
}
print(count-1)

g = upset(upset_mtx_pan, intersect = colnames(upset_mtx_pan), min_size=10,
          sort_intersections_by = 'degree', sort_intersections = 'ascending',
          mode = 'inclusive_intersection', width_ratio = 0.1, 
          min_degree=2, set_sizes = FALSE,
          queries = uquery,
          base_annotations=list(
            'Intersection size'=intersection_size(
              mode = 'inclusive_intersection',
              text=list(vjust=-0.2,hjust=0.5,angle=0,size=3),
              #mapping=aes(fill=mpaa),
              text_colors=c(on_background='black', on_bar='white'))
            + annotate(
              geom='text', x=Inf, y=Inf, color='#FF7F0EFF',
              label='Significant overlap', size=5,
              vjust=1, hjust=1
            ))
          # set_sizes = upset_set_size() + 
          #   scale_y_continuous(breaks=c(4000,2000,0))
) & theme(plot.background=element_rect(fill='transparent', color=NA))
g
ggsave(file.path(plot_out, 'panCan_pos_upset.pdf'),
       plot=g, width=11.69*1,height=8.27*0.5,units='in',device='pdf',dpi=300)

### Get top five GO for each significant overlap ====
top5terms = list()
ids = c()
count = 0
for (i in uquery){
  its = i$intersect
  ids = c(ids, its)
  id = paste(its, collapse =' X ')
  genes = rownames(upset_mtx_pan)[which(rowSums(upset_mtx_pan[,its]) == length(its))]
  if (length(genes)>=20){
    count = count + 1
    print(id)
    test = do_GO(genes, background = hs$gene[which(hs$beta>0)])
    ps = pcs_GO_out(test, filename = paste('panCan_posOvl_GO',id,'.pdf',sep=''), 
                    dir = plot_out, cneplt = F,
                    size = c(8.27*0.7,11.69*0.5))
    
    rs = test@result
    rs = rs[rs$p.adjust < 0.05,]
    if (nrow(rs)>0){
      top5terms[[id]] = test@result$Description[1:min(5,nrow(test@result))]
    } else {
      top5terms[[id]] = NULL
    }
  } else {
    print('Fewer than 20')
    print(id)
  }
}
top5term_df = as.data.frame(sort(table(unlist(top5terms)), decreasing = F))
top5term_df$Var1 = factor(top5term_df$Var1, levels = top5term_df$Var1)
g = ggplot(top5term_df[(nrow(top5term_df)-5):nrow(top5term_df),]) + scale_x_continuous(expand = c(0,0)) +
  geom_bar(aes(y=Var1,x=Freq/count), stat = 'identity') + mytme + 
  labs(x='Frequency', y='') +
  geom_text(aes(label=Freq, y=Var1, x=Freq/count), nudge_x = -0.1, color='white') +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0), breaks=seq(0,1,0.5))
ggsave(file.path(plot_out, 'panCan_pos_upset_GO.pdf'),
       plot=g, height=11.69*0.3,width=8.27*0.5,units='in',device='pdf',dpi=300)



## print overlap among cancers ====
### positive genes ====
df_ovl = data.frame()
for (i in 1:ncol(upset_mtx_pan)){
  tar = rownames(upset_mtx_pan[rowSums(upset_mtx_pan) >= i,])
  df_ovl = rbind(df_ovl, c(i,length(tar),length(intersect(tar, pos.ctrs))))
}
colnames(df_ovl) = c('trh', 'tarGene', 'ctrGene')
df_ovl = gather(df_ovl, key='gp', value='Count', 2:3)
df_ovl$gp = factor(df_ovl$gp, levels = c('tarGene', 'ctrGene'))
g = ggplot(df_ovl, aes(x=trh,y=Count, color=gp)) +
  geom_vline(xintercept = 3, color='red', size=1) +
  geom_line(size=1.5) + geom_point(size=1.5, color='black') + mytme +
  scale_color_d3(palette = 'category20', name='', labels = c('Identified positive genes', 'Positive control genes')) +
  labs(x=str_wrap('Identified in more than # cancer types', width = 20)) +
  scale_y_log10() +
  theme(legend.position = c(0.75,0.9)) +
  scale_x_continuous(breaks = 1:ncol(upset_mtx_pan))
ggsave(file.path(plot_out, 'panCan_pos_trh.pdf'),
       plot=g, height=11.69*0.4,width=8.27*0.8,units='in',device='pdf',dpi=300)
### negative genes ====
df_ovl = data.frame()
for (i in 1:ncol(upset_mtx_pan)){
  tar = rownames(upset_mtx_pan[rowSums(upset_mtx_pan) >= i,])
  df_ovl = rbind(df_ovl, c(i,length(tar),length(intersect(tar, neg.1.ctrs)),
                           length(intersect(tar, neg.2.ctrs))))
}
colnames(df_ovl) = c('trh', 'tarGene', 'ctrGene.1', 'ctrGene.2')
df_ovl = gather(df_ovl, key='gp', value='Count', 2:4)
df_ovl$gp = factor(df_ovl$gp, levels = c('tarGene', 'ctrGene.1', 'ctrGene.2'))
g = ggplot(df_ovl, aes(x=trh,y=Count, color=gp)) +
  geom_vline(xintercept = 3, color='red', size=1) +
  geom_line(size=1.5) + geom_point(size=1.5, color='black') + mytme +
  scale_color_d3(palette = 'category20', name='', 
                 labels = c('Identified negative genes', 'Negative control genes 1', 'Negative control genes 2')) +
  labs(x=str_wrap('Identified in more than # cancer types', width = 20)) +
  scale_y_log10() +
  scale_x_continuous(breaks = 1:ncol(upset_mtx_pan)) + 
  theme(legend.position = c(0.75,0.9))
ggsave(file.path(plot_out, 'panCan_neg_trh.pdf'),
       plot=g, height=11.69*0.4,width=8.27*0.8,units='in',device='pdf',dpi=300)

### Run GO on degree x genes ====
deg_trh = 3
tar = rownames(upset_mtx_pan[rowSums(upset_mtx_pan) >= deg_trh,])
# !!! change sign of the background !!!
test = do_GO(tar, background = unique(hs$gene[which(hs$beta>0)])) # tar: 234 genes degree 3, 74 genes degree 4
ps = pcs_GO_out(test, filename = 'panCan_pos_trh3_GO.pdf', dir = plot_out, cneplt = F,
                size = c(8.27*0.7,11.69*0.5)) # width, height
test_rvg = run_revigo(test)
test_rvg = test[test@result$ID%in%test_rvg$`Term ID`,]
repair_genes = strsplit(test_rvg['GO:0000278','geneID'],split='/')[[1]]
repair_ctr = load_controls(upset_mtx_pan[repair_genes,])
repair_ctr = repair_ctr[repair_genes,]

mean_beta = hs[hs$gene %in% tar,] %>% 
  group_by(gene) %>% summarise(mbeta=mean(beta))

test = do_GO(tar, background = unique(hs$gene[which(hs$beta<0)])) # tar: 146 genes, 49 genes degree 4
ps = pcs_GO_out(test, filename = 'panCan_neg_trh3_GO.pdf', dir = plot_out, cneplt = F,
                size = c(8.27*0.7,11.69*0.5)) # p53 markers


# Test hot spot genes in BRCA in CCLE ========
tcga.brca = subset(df_coll, df_coll$protein_change=='hot_spot' & df_coll$experiment=='tcga_brca_raw_seq')
tcga.up = tcga.brca$gene[tcga.brca$beta>0]

## Load controls ====
mylist = load_controls(ng_ctrs = T)
names(mylist) = c('CCLE RNA up', 'CCLE RNAi down', 'CCLE RNA NS up','CCLE CRISPR down', 
                  'known_wt_up', 'known_wt_down', 'R273H ChIP peaks',
                  'ng1', 'ng2', 'Positive controls')
mylist[['TCGA RNA up']] = tcga.up
mylist = mylist[setdiff(names(mylist), c('ng1', 'ng2'))]

um = gen_upSet_mtx(mylist)
um = um[um$known_wt_up==0 & um$known_wt_down==0,]
um = um[,-which(colnames(um) %in% c('known_wt_up', 'known_wt_down'))]
um = um[rowSums(um)>=1,]
um_plt = um[um$`TCGA RNA up`>0,] # filter identified in eQTL only
print(colSums(um))

### Shuffle upset matrix, how many observed is significant? ====
# T: must be 1, F: must be 0, NA: no matter. 
test_condition = c(NA, T, T, NA, NA, NA, T)
ovl_p = run_overlap_test(um, test_condition = test_condition, n_loop = 1000)

### Plot upset ====
# https://krassowski.github.io/complex-upset/articles/Examples_R.html
g = upset(um_plt, colnames(um_plt), min_size = 1, width_ratio = 0.15,
      min_degree = 2,
      sort_intersections_by = c('degree'), sort_intersections = c('ascending'),
      queries = list(upset_query(set='Positive controls', color='#2CA02CFF', fill='#2CA02CFF'),
                     upset_query(set='TCGA RNA up', color='#1F77B4FF', fill='#1F77B4FF'),
                     upset_query(intersect=c('TCGA RNA up', 'R273H ChIP peaks'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'CCLE RNA NS up'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'CCLE RNA up'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'CCLE RNA up', 'CCLE RNA NS up'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'R273H ChIP peaks', 'CCLE RNAi down'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'R273H ChIP peaks', 'CCLE RNA up', 'CCLE RNA NS up'), fill='#FF7F0EFF', color='#FF7F0EFF'),
                     upset_query(intersect=c('TCGA RNA up', 'R273H ChIP peaks', 'CCLE RNA up'), fill='#FF7F0EFF', color='#FF7F0EFF')),
      base_annotations=list(
        'Intersection size'=intersection_size(
          mode = 'inclusive_intersection',
          text=list(vjust=-0.2,hjust=0.5,angle=0,size=3),
          text_colors=c(on_background='black', on_bar='white'))
          + annotate(
            geom='text', x=Inf, y=Inf, color='#FF7F0EFF',size=5,
            label='Intersections with positive controls',
            vjust=1, hjust=1
          ))
      )
g
ggsave(file.path(plot_out, 'BRCA', 'upset.pdf'),
       plot=g, width=11.69*1,height=8.27*0.5,units='in',device='pdf',dpi=300)


### Enrichment analysis ====
# load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
# mtx = as.data.frame(mtx)
# bg = rownames(mtx)[which(mtx$BRCA==1)]
bg = tcga.brca$gene[tcga.brca$beta > 0] # use all identified as bg
print(colnames(um))
rownames(tcga.brca) = tcga.brca$gene
tcga.brca = tcga.brca[order(tcga.brca$beta, decreasing = T),]
# ChIP and TCGA
gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,NA,NA,NA,T,NA,T))] # 584 genes

## gsea
to_gsea = tcga.brca[gene,'beta']
names(to_gsea) = gene
to_gsea = sort(to_gsea, decreasing = T)
# plot(to_gsea)
test_kegg = do_GSEA(to_gsea, pvalue=1)
rs = test_kegg@result[which(test_kegg@result$qvalues<0.25),]
gseaplot2(test_kegg,1)

test = do_GO(gene, background = bg, ont='BP')
sign1 = ext_gene_GO(test@result$geneID[1:5], do_intersect = F)
g = dotplot(test) + labs(title= 'TCGA RNA up X ChIP-seq peaks')
ggsave(file.path(plot_out,'BRCA', '1_chip-tcga_GO.pdf'),
       width=6, height=4, units='in', device='pdf', dpi=300, plot = g)
# TCGA RNA and CCLE RNA
gene = rownames(um)[get_idx_condt(um, test_condition = c(T,NA,NA,NA,NA,NA,T))] # 565 genes
test = do_GO(gene, background = bg, ont='BP')
g = dotplot(test) + labs(title= 'TCGA RNA up X CCLE RNA up') # nothing
# TCGA RNA and CCLE RNA NS
gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,NA,T,NA,NA,NA,T))] # 177 genes
test = do_GO(gene, background = bg, ont='BP')
sign3 = ext_gene_GO(test@result$geneID[1:5], do_intersect = F)
g = dotplot(test) + labs(title= 'TCGA RNA up X CCLE RNA NS up')
ggsave(file.path(plot_out,'BRCA', '3_ccleNS-tcga_GO.pdf'),
       width=6, height=4, units='in', device='pdf', dpi=300, plot = g)
# ChIP and TCGA, RNAi
gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,T,NA,NA,T,NA,T))] # 13 genes
test = do_GO(gene, background = bg, ont='BP')
g = dotplot(test) + labs(title= 'TCGA RNA up X ChIP-seq peaks X CCLE RNAi down') # nothing
# ChIP and TCGA, CCLE RNA
gene = rownames(um)[get_idx_condt(um, test_condition = c(T,NA,NA,NA,T,NA,T))] # 109 genes
test = do_GO(gene, background = bg, ont='BP')
wnt_genes = test@result$geneID[1]
g = dotplot(test) + labs(title= str_wrap('TCGA RNA up X ChIP-seq peaks X CCLE RNA up', width=28))
ggsave(file.path(plot_out,'BRCA', '5_chip-tcga-CCLE_GO.pdf'),
       width=6, height=4, units='in', device='pdf', dpi=300, plot = g)
# TCGA and CCLE RNA and CCLE RNA ns
gene = rownames(um)[get_idx_condt(um, test_condition = c(T,NA,T,NA,NA,NA,T))] # 91 genes
test = do_GO(gene, background = bg, ont='BP') # nothing
g = dotplot(test) + labs(title= str_wrap('TCGA RNA up X CCLE RNA up X CCLE NS RNA up', width = 28))
# ChIP and TCGA, CCLE RNA, CCLE RNA NS
gene = rownames(um)[get_idx_condt(um, test_condition = c(T,NA,T,NA,T,NA,T))] # 26 genes
test = do_GO(gene, background = bg, ont='BP') # nothing
g = dotplot(test) + labs(title= str_wrap('TCGA RNA up X ChIP-seq peaks X CCLE RNA up X CCLE RNA NS up', width=28))

### Profile microtubule genes ====
sign3_mtx = um[rownames(um) %in% sign3,]
rownames(sign3_mtx[sign3_mtx$`CCLE RNA up`==1,])
sign1_mtx = um[rownames(um) %in% sign1,]
rownames(sign1_mtx[sign1_mtx$`CCLE RNA up`==1 & sign1_mtx$`CCLE RNA NS up`==1,])

### Profile wnt signalling genes ====
# wnt_genes = strsplit(wnt_genes, split='/')[[1]]
## from wnt signature (old)
GOI = c("FZD9","GPRC5B","RUVBL1","CSNK2A1","RPS12","RNF220","TNFAIP3")
## top 5 rna regulation terms (sign1)
GOI = c("BYSL","NCL","PPIH","THOC5")
## top 5 microtubule & mitosis terms (sign3)
GOI = c("ARHGEF2","CDC20","CHEK2","CRYAB","KIF4A","KIFC3","MAP7D3","MYH9","RCC1", "STIL","SUN2")
## sign1 and sign3 with ChIP peak, CCLE NS
GOI = c("BYSL","NCL","PPIH","THOC5",'ARHGEF2', 'CDC20')

dt.tcga = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData', 
                          ann_bin_mut_list = c('hot_spot'), mode='tcga')
dt.ccle = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd/BRCA/clean_data.RData',
                          ann_bin_mut_list = c('hot_spot'), mode = 'ccle')
df.tcga = cbind(t(assay(dt.tcga[[1]][GOI,,'RNA'])), as.data.frame(dt.tcga[[1]]@colData[,c('p53_state', 'has_hot_spot')]))
df.ccle = cbind(t(assay(dt.ccle[[1]][GOI,,'RNA'])), as.data.frame(dt.ccle[[1]]@colData[,c('p53_state', 'has_hot_spot')]))
# !!! choose one !!!
df.plt = df.ccle
df.plt = df.tcga
df.plt = gather(df.plt %>% mutate(sample=rownames(df.plt)), key='gene', value='expression', 1:length(GOI))
df.plt$gene = factor(df.plt$gene, levels=GOI)
df.plt = df.plt[df.plt$p53_state == 'Wildtype' | (df.plt$p53_state=='missense' & df.plt$has_hot_spot==1) |
                  (df.plt$p53_state == 'nonsense'),]
df.plt = df.plt[order(df.plt$gene, df.plt$has_hot_spot),]
df.plt$sample = factor(df.plt$sample, levels = df.plt$sample[1:(nrow(df.plt)/length(GOI))])
myPalette = colorRampPalette(c("#1F77B4FF",'white', '#D62728FF'))
sc = scale_fill_gradientn(colours = myPalette(50), name='Normalized expression')
tb = c(sum(df.plt$p53_state=='Wildtype'), sum(df.plt$has_hot_spot==1), sum(df.plt$p53_state=='nonsense'))
df.plt$has_hot_spot[df.plt$has_hot_spot==0] = paste('WT', tb[1] / length(GOI), sep=' n=')
df.plt$has_hot_spot[df.plt$has_hot_spot==1] = paste('HS', tb[2] / length(GOI), sep=' n=')
df.plt$has_hot_spot[df.plt$p53_state=='nonsense'] = paste('NS', tb[3]/length(GOI), sep=' n=')
g = ggplot(df.plt) +
  geom_tile(aes(x=sample, y=gene, fill=expression), width=1) +
  #geom_vline(xintercept = n_sample[1], size=1, alpha=0.7) +
  sc + mytme + labs(x='', y='') +
  facet_grid(~has_hot_spot,scale='free_x', space='free') +
  theme(legend.direction = 'horizontal', axis.text.x = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold', size=12),
        panel.border = element_rect(color='black', fill='transparent', size=1.5)) 
# CCLE width: 5, TCGA width: 10
ggsave(file.path(plot_out,'BRCA', 'WntHeatmapTCGA.pdf'),
       plot = g + theme(legend.key.width = unit(0.4,units = 'inch')),
       width=8.27*1.5, height=4, units='in', device='pdf', dpi=300, bg = 'transparent')
ggsave(file.path(plot_out,'BRCA', 'WntHeatmapCCLE.pdf'),
       plot = g + theme(legend.key.width = unit(0.3,units = 'inch')),
       width=8.27*0.6, height=4, units='in', device='pdf', dpi=300, bg = 'transparent')
