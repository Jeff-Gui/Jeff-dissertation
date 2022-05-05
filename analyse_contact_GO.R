## There is the neurogenesis and cell adhesion signature in eQTL enriched for contact mutants,
## Are they expressed in the neural system?
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult TCGA-pan_VS-wt
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'coreVScontact')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../set_theme.R')
source('../revigo_utils.R')
library(tidyverse)
library(stringr)
library(ggvenn)
library(enrichplot)
library(ggsci)
library(ComplexUpset)
source('../overlap_utils.R')
exp_plt_out = file.path(dir_home, 'plots', 'coreVScontact')

# Load data ====
## load eQTL
div = read.table(file.path(data_out, 'BRCA_division_qtl.txt'), header = T)
mig = read.table(file.path(data_out, 'BRCA_migration_qtl.txt'), header = T)
nev = read.table(file.path(data_out, 'BRCA_neural_qtl.txt'), header = T)
brca = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData',
                       ann_bin_mut_list = c('hot_spot', 'sandwich', 
                                            'contact', 'conformation'))
## load Xena
load('/Users/jefft/Desktop/p53_project/datasets/Xena/BRCA_GTEX_TCGA/TCGA-BRCA_GTEX-nerve-brain-breast_NormCount.RData')
phen = read.table('/Users/jefft/Desktop/p53_project/datasets/Xena/BRCA_GTEX_TCGA/TCGA-BRCA_GTEX-nerve-brain-breast_phenotype.txt', sep='\t', header = T)

# Extract gene expr ====
GOI = unique(div$gene)
gtex = phen$sample[phen$X_study=='GTEX']
co_sample = intersect(rownames(brca[[1]]@colData), colnames(toil_sel))
brca[[1]] = subsetByColumn(brca[[1]], co_sample)
table(brca[[1]]@colData$p53_state)
expr = toil_sel[GOI,c(gtex, rownames(brca[[1]]@colData))]
expr = na.omit(expr) # 5 genes missing
print(paste('Input gene:', length(GOI)))
print(paste('Overlapping dim:',paste(dim(expr), collapse = '-')))
GOI = rownames(expr)

# Distribution ====
## there are samples with really low expression in the nev gene
hc = hclust(dist(t(expr)))
tree = cutree(hc, k=2)
table(tree) # 5 GTEX samples
bad_sample = names(tree[which(tree==2)]) # remove them from analysis
expr = expr[,!colnames(expr) %in% bad_sample]
ann = phen[phen$sample %in% colnames(expr),c('sample','X_primary_site','X_study')]
rownames(ann) = ann$sample
ann = ann[,-1]
ann$X_primary_site[ann$X_primary_site=='Breast' & ann$X_study=='TCGA'] = 'Breast cancer'
ann$p53_state = 'wt'
ann[rownames(brca[[1]]@colData)[brca[[1]]@colData$p53_state!='Wildtype'],'p53_state'] = 'mut'
ann[rownames(brca[[1]]@colData)[brca[[1]]@colData$has_contact==1],'p53_state'] = 'contact'
ann = ann[order(ann$p53_state),]
ann$subtype=NA
ann[rownames(brca[[1]]@colData), 'subtype'] = brca[[1]]@colData$SUBTYPE
ann$subtype[ann$subtype==''] = NA
table(ann$p53_state)
table(ann$subtype)

expr_z_score = as.data.frame(t(scale(t(expr))))
expr_z_score = expr_z_score[,rownames(ann)]
pheatmap::pheatmap(expr[,rownames(ann)[ann$X_study=='TCGA' | ann$X_primary_site %in% c('Breast')]], 
                   show_rownames = F, show_colnames = F, 
                   annotation_col = ann[,c('X_primary_site', 'p53_state', 'subtype')],
                   cluster_rows = T, cluster_cols = F, scale = 'row',
                   color = colorRampPalette(c('blue','white','red'))(100),
                   filename = file.path(data_out, 'division_heatmap.pdf'))

# PCA ====

# Compare cancer to normal for gene in eQTL ====
qtl = read.table(file.path(eqtl_out, 'tcga_brca_raw_seq', 'trans_eqtl_fdr005.txt'), sep='\t', header = T)
qtl = subset(qtl, qtl$beta > 0 & qtl$protein_change %in% c('sandwich', 'conformation', 'contact'))

# Run DEG ====
## Plain method to generate DEG ====
ps = c()
fc = c()
# gene = intersect(unique(qtl$gene), rownames(toil_sel))
normal_breast = rownames(ann)[ann$X_study=='GTEX' & ann$X_primary_site=='Breast']
cancer_breast = rownames(ann)[ann$X_study=='TCGA']
for (i in 1:nrow(toil_sel)){
  a = as.numeric(toil_sel[i,normal_breast])
  b = as.numeric(toil_sel[i,cancer_breast])
  fc = c(fc, mean(b) - mean(a))
  ps = c(ps, wilcox.test(a,b)$p.value)
  if (i %% 5000==0){
    print(i)
  }
}
deg = data.frame('gene'=rownames(toil_sel), 'foldchange'=fc, 'pvalue'=ps, 'padj'=p.adjust(ps, method = 'fdr'))
write.table(deg, file.path(data_out, 'TCGA-BRCA_p53mut_VS_XenaToil_Breast_DEG.txt'), sep='\t', row.names=F, quote = F)

## DEG from eQTL (missense p53 TCGA VS breast normal GTEX) ====
deg = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/test/outputs/tcga_gtex_brca/trans_eqtl_fdr005.txt', sep='\t', header = T)
colnames(deg) = c('snps','gene','stat','pvalue','padj', 'foldchange','protein_change')
deg = deg[deg$protein_change=='all',]
deg = deg[,c('gene','pvalue','padj','foldchange')]
if (anyDuplicated(deg$gene)){
  deg = deg[-which(duplicated(deg$gene)),]
}
deg_test = read.table(file.path(data_out, 'TCGA-BRCA_p53mut_VS_XenaToil_Breast_DEG.txt'), sep='\t', header = T)
miss = setdiff(deg_test$gene, deg$gene)
aug = as.data.frame(matrix(rep(c(NA,1,1,0), length(miss)),nrow=length(miss), byrow = T))
colnames(aug) = colnames(deg)
aug$gene = miss
deg = rbind(deg, aug)


## Summarise deg ====
rownames(deg) = deg$gene
qtl_nm = left_join(qtl, deg, by='gene')
core = as.data.frame(table(qtl_nm$gene))
unq = core[core$Freq==1,]$Var1
core = core[core$Freq==3,]$Var1
rest = unique(qtl_nm$gene)
rest = rest[!rest %in% core]
contact = intersect(unique(qtl_nm[qtl_nm$protein_change=='contact', 'gene']), unq)
sandwich = intersect(unique(qtl_nm[qtl_nm$protein_change=='sandwich', 'gene']), unq)
conform = intersect(unique(qtl_nm[qtl_nm$protein_change=='conformation', 'gene']), unq)

core_go = do_GO(core, background = unique(qtl$gene))
dotplot(core_go)
core_div = ext_gene_GO(core_go@result$geneID[1:5])

annotate = function(x){
  x = as.character(x)
  # x: a vector of genes
  deg2 = deg
  missing = x[!x%in%deg$gene]
  app = as.data.frame(matrix(NA, nrow=length(missing), ncol=ncol(deg)))
  colnames(app) = colnames(deg)
  deg2 = rbind(deg2, app)
  rownames(deg2) = c(rownames(deg),missing)
  sub = deg2[x,]
  sub$ann = 'ND'
  sub$ann[sub$foldchange<0] = 'down'
  sub$ann[sub$foldchange>0] = 'up'
  sub$ann[sub$padj>=0.05] = 'noSig'
  return(sub[x,'ann'])
}

smm_list = list(core, core_div, rest, contact, unique(nev$gene), unique(mig$gene),
                sandwich, conform, unique(qtl$gene))
nsize = sapply(smm_list, length)
smm = lapply(smm_list,
       function(x){return(table(annotate(x))/length(x))})
smm = Reduce(cbind, smm)
nm = c('core','div','nonCore', 'contact','nev','mig', 'sandwich', 'conformation', 'union')
flv = c('union','core','div', 'nonCore', 'contact','nev' ,'mig','conformation', 'sandwich')
colnames(smm) = nm
names(nsize) = nm
smm = gather(as.data.frame(smm) %>% mutate(cls=rownames(smm)), key='gp', value='freq',1:ncol(smm)) %>%
  mutate(gp=factor(gp, levels = flv))
## shuffle and test significance of down
## - randomly choose same number of genes from all eQTL input, get from DEG table
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
qtl_pool = rownames(mtx)[mtx$BRCA==1]

library(parallel)
#detectCores()
cl = makeCluster(8)
clusterExport(cl, c('deg','qtl_pool','nsize','annotate'))

n_loop = 1000
null_dis = parLapply(cl, 1:n_loop, function(x){
  ys = list()
  for (i in 1:length(nsize)){
    gn = sample(qtl_pool, nsize[i], replace = F)
    atn = annotate(gn)
    y = table(atn) / length(gn)
    y = y[c('down','ND','noSig','up')]
    y['gp'] = names(nsize[i])
    ys[[i]] = y
  }
  return(Reduce(rbind, ys))
})
null_dis = Reduce(rbind, null_dis)
stopCluster(cl)
colnames(null_dis) = c('down','ND','noSig','up', 'gp')
rownames(null_dis) = NULL
null_dis = as.data.frame(null_dis)
null_dis = null_dis[order(null_dis$gp),]
print(sum(is.na(null_dis)))
null_dis[is.na(null_dis)] = 0

lb_gp = c('Union', 'Core','Core (cell cycle GO)', 'Non-core', 'Contact','Contact (neural GO)','Contact (adhesion GO)', 'Conformation', 'Sandwich')
g = ggplot(null_dis %>% gather(key='cls', value='freq', 1:4) %>%
         mutate(gp=factor(gp, levels = rev(flv)))) +
  geom_boxplot(aes(y=gp,x=100*as.numeric(freq)),alpha=1, width=0.5, outlier.size = 0.5) +
  geom_point(data = smm, aes(x=100*freq,y=gp), color='red', shape=15) +
  facet_wrap(~cls, scales = 'free_x', nrow=1,
             labeller = labeller(cls=c('down'='Low','ND'='Not detected', 
                                       'noSig'='No significance', 'up'='High'))) +
  scale_y_discrete(labels=rev(lb_gp)) +
  scale_x_continuous(expand = c(0.1,0.1)) +
  labs(x='Frequency (%)', y='') +
  mytme +
  theme(legend.direction = 'horizontal', strip.background = element_rect(fill='transparent'),
        strip.text = element_text(size=12, face='bold'),
        panel.border = element_rect(fill='transparent',size=1)) +
  scale_color_d3(palette = 'category20')

# library(gridExtra)
# p_tab = tableGrob(paste('n=',nsize[rev(flv)],sep=''))
# grid.arrange(g, p_tab, ncol = 2, padding=0)

ggsave(file.path(plot_out,'BRCA', 'compositionVSnormal_bootstrap_fromQTL.pdf'),
       plot = g, width=11.67*0.8, height=8.27*0.5, units='in', device='pdf', dpi=300, bg = 'transparent')


g = ggplot(smm, 
       aes(x=gp, y=freq, fill=cls)) +
  geom_bar(stat = 'identity') + mytme +
  scale_fill_d3(palette = 'category20', name='', labels=c('Low', 'Not detected', 'No significance', 'High')) +
  labs(x='', y='Frequency') +
  geom_hline(yintercept = 0.5, color='black', linetype='dotted', size=1) +
  scale_x_discrete(labels=lb_gp) +
  geom_text(data = data.frame('gp'=names(nsize),'size'=nsize), aes(x=gp,y=0.07,label=paste('n=',size,sep=''),fill=NA),
     size=3, color='white',fontface='italic') +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.direction = 'horizontal', axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = 'top')
ggsave(file.path(plot_out,'BRCA', 'compositionVSnormal_fromQTL.pdf'),
       plot = g, width=8.27*0.8, height=5, units='in', device='pdf', dpi=300, bg = 'transparent')

# plot beta in two eQTL comparisons ====
library(ggrepel)
lb_gene = unique(mig[mig$peak.over.chek1==1 & mig$beta>0,]$gene)

two_qtl = na.omit(qtl_nm)
# two_qtl = two_qtl[two_qtl$padj<0.05,]
two_qtl$label = two_qtl$gene
two_qtl$label[!two_qtl$label %in% lb_gene] = NA
# filter for contact-specific genes
two_qtl_plt = two_qtl[two_qtl$protein_change=='contact' & two_qtl$gene %in% contact,]
g = ggplot(two_qtl_plt, aes(y=beta, x=foldchange)) +
  geom_point(size=0.1) +
  # geom_bin2d() +
  # facet_wrap(~protein_change, ncol=1) +
  geom_point(data = two_qtl_plt[!is.na(two_qtl_plt$label),], color='red') +
  geom_text_repel(aes(label=label), box.padding = 0.1) +
  labs(y='beta (p53 contact Mut VS p53 WT)', x='beta (all p53 missense Mut VS normal tissue)') +
  mytme
ggsave(file.path(plot_out,'BRCA', 'contactSpc_beta.pdf'),
       plot = g, width=8.27*0.7, height=11.67*0.4, units='in', device='pdf', dpi=300, bg = 'transparent')



# plot per gene (Toil UQ original data) ====
sample_sel = rownames(ann)[ann$X_study=='TCGA' | ann$X_primary_site %in% c('Breast')]
gene = 'PPM1F'
df.plt = cbind(t(toil_sel[gene,sample_sel]), ann[sample_sel,])
df.plt$itg_gp = paste(df.plt$X_primary_site, df.plt$p53_state, sep='_')
# df.plt[[gene]] = (df.plt[[gene]] - mean(df.plt[[gene]])) / sd(df.plt[[gene]])
ggplot(df.plt, aes(x=itg_gp, y=get(gene))) +
  geom_violin() + geom_boxplot(width=0.2, outlier.shape = NA) +
  labs(x='', y='Toil UQ', title=gene) +
  stat_compare_means(comparisons = list(c('Breast cancer_contact', 'Breast cancer_wt'),
                                        c('Breast cancer_wt', 'Breast_wt'))) +
  mytme

# plot per gene (rankNorm data - eQTL input) ====
fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA-GTEX/clean_data.RData'
tcga.gtex = load_clean_data(fp,ann_bin_mut_list = c('contact'))
meta = as.data.frame(tcga.gtex[[1]]@colData)
meta$p53_state[meta$has_contact==1] = 'contact'
gene = 'SUN2'
df.plt = cbind(t(assay(tcga.gtex[[1]][gene,rownames(meta),'RNA'])), meta)
df.plt$itg_gp = paste(df.plt$study, df.plt$p53_state, sep=' ')
df.plt$itg_gp = recode_factor(df.plt$itg_gp,'GTEX Wildtype'='nm',
                              'TCGA missense' = 'cmis',
                              'TCGA contact' = 'ccon',
                              .ordered=TRUE)
ggplot(df.plt, aes(x=itg_gp, y=get(gene))) +
  geom_violin() + geom_boxplot(width=0.2, outlier.shape = NA) +
  labs(x='', y='Normalised expression', title=gene) +
  scale_x_discrete(labels = c('Normal breast',
                              'BRCA\np53 non-contact\nmissense',
                              'BRCA\np53 contact\nmissense')) +
  stat_compare_means(comparisons = list(c('nm', 'cmis'),
                                        c('cmis', 'ccon'),
                                        c('nm', 'ccon'))) +
  # stat_compare_means(method = 'anova') +
  mytme

# plot heatmap ====
df.plt = as.data.frame(t(expr))
df.plt = gather(df.plt %>% mutate(sample=rownames(df.plt)), key='gene', value='expression', 1:length(GOI))
df.plt$gene = factor(df.plt$gene, levels=GOI)
idt = paste('has_', mut_groups, sep='')
all_mis = intersect(which(df.plt$p53_state=='missense'), which(rowSums(df.plt[idt])==1))
df.plt = df.plt[c(which(df.plt$p53_state == 'Wildtype'),
                  which(df.plt$p53_state == 'nonsense'),
                  all_mis),]
df.plt = df.plt[order(df.plt$gene, rowSums(df.plt[idt])),]
df.plt$sample = factor(df.plt$sample, levels = df.plt$sample[1:(nrow(df.plt)/length(GOI))])
myPalette = colorRampPalette(c("#1F77B4FF",'white', '#D62728FF'))
sc = scale_fill_gradientn(colours = myPalette(50), name='Normalized expression')

df.plt$plt_gp = df.plt$p53_state
df.plt$plt_gp = gsub('Wildtype', paste('WT', sum(df.plt$p53_state=='Wildtype') / length(GOI), sep=' n='), df.plt$plt_gp)
df.plt$plt_gp = gsub('nonsense', paste('NS', sum(df.plt$p53_state=='nonsense') / length(GOI), sep=' n='), df.plt$plt_gp)
for (i in idt){
  lb = paste(str_to_sentence(strsplit(i, split='_')[[1]][2]),sum(df.plt[i]==1) / length(GOI), sep=' n=')
  df.plt$plt_gp[df.plt[i]==1] = lb
}

g = ggplot(df.plt) +
  geom_tile(aes(x=sample, y=gene, fill=expression), width=1) +
  #geom_vline(xintercept = n_sample[1], size=1, alpha=0.7) +
  sc + mytme + labs(x='', y='') +
  facet_grid(~plt_gp,scale='free_x', space='free') +
  theme(legend.direction = 'horizontal', axis.text.x = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold', size=12),
        panel.border = element_rect(color='black', fill='transparent', size=1.5)) 

