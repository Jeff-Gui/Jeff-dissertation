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
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)
library(ggvenn)
library(enrichplot)
library(ComplexUpset)

gl = load_GO_out(exp_home = dir_home, gene_ratio_min = 0)
go_df = gl$df
go_df = go_df[grep('GO', go_df$ID),]

# Core VS Contact mutant, GO =====================================
signn = 'neg'
exp_plt_out = file.path(dir_home, 'plots', 'coreVScontact')

cc_go_df = subset(go_df, go_df$mutation %in% c('conformation', 'contact', 'sandwich'))
cc_go_df = cc_go_df[grep('GO', cc_go_df$ID),]

cc_go_df = subset(cc_go_df, cc_go_df$sign==signn)
table(cc_go_df$experiment)

plt.list = list()
for (i in unique(cc_go_df$cancer)){
  test = cc_go_df[cc_go_df$cancer==i,]
  g = ggvenn(list('conformation'=test$ID[test$mutation=='conformation'], 
                  'contact'=test$ID[test$mutation=='contact'],
                  'sandwich'=test$ID[test$mutation=='sandwich'])) +
    labs(title=i)
  plt.list[[i]] = g
}
plt.list %>% marrangeGrob(ncol=3, nrow=3, top = '') %>%
  ggsave(file.path(exp_plt_out, 'GO_overview_neg.pdf'),
         plot=., width=10,height=10,units='in',device='pdf',dpi=300)

## GO enrichment for core mutation found in most cancers ==========
signn = 'neg'
core_gos = go_df[go_df$mutation == 'conformation',]
core_gos_pos = core_gos[core_gos$sign==signn,]
coll = list()
for (i in unique(core_gos_pos$cancer)){
  coll[[i]] = core_gos_pos$ID[which(core_gos_pos$cancer == i)]
}
g = ggvenn(coll, stroke_color = 'white')
g

#### dot plot of each top 5 terms for each group
top_go_plt = data.frame()
for (i in unique(core_gos_pos$cancer)){
  sub_go_df = subset(core_gos_pos, core_gos_pos$cancer==i & core_gos_pos$sign == signn)
  if (nrow(sub_go_df)>0){
    top_go_plt = rbind(top_go_plt, sub_go_df[order(sub_go_df$p.adjust)[1:(min(10,nrow(sub_go_df)))],])
  }
}

# top_go_plt$Description = factor(top_go_plt$Description, levels = unique(top_go_plt$Description))
g = ggplot(top_go_plt, aes(x=gene_ratio_val)) +
  geom_point(aes(y=Description, color=-log10(p.adjust), size=Count, group=cancer)) +
  facet_wrap(~cancer, scale='free', nrow = 2) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=25)) +
  theme(axis.text.x = element_text(angle=0, vjust = 1, size=16),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size=16, face='bold')) +
  mytme
g
g %>% ggsave(file.path(exp_plt_out, paste(signn, '_pan_conformation_GO_dotPlot.pdf', sep='')),
             plot=., width=25,height=15,units='in',device='pdf',dpi=300)


# There is only GO enrichment of contact mutant in BRCA and LGG.
# load CCLE
load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
ccle.core = test_ccle(dt, check_mut = 'core', genes=NULL, check_site = 'Breast')
ccle.contact = test_ccle(dt, check_mut = 'contact', genes=NULL, check_site = 'Breast')
ccle.core = ccle.core$ccle.valid
ccle.contact = ccle.contact$ccle.valid
save(ccle.core, ccle.contact, file = file.path(data_out, 'ccle_breast_ContactCorevsWT_allgenes.RData'))

# Gene ================
## Load dt ====
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
# load(file.path(data_out, 'ccle_breast_ContactCorevsWT_allgenes.RData'))
beta_cutoff=0
coll = load_eQTL_output(eqtl_out, beta=beta_cutoff, exclude = 'tcga_nine_pool')
df_coll = data.frame()
for (i in names(coll)){
  nm = toupper(strsplit(i, split='_')[[1]][2])
  df = coll[[i]]
  if (class(df)!='logical'){
    df = subset(df, abs(df$beta) > beta_cutoff)
    df_coll = rbind(df_coll, df)
  }
}
df_coll = subset(df_coll, df_coll$protein_change %in% c('contact', 'conformation', 'sandwich'))

## Overlap overview in cancers ====
plt.list = list()
for (i in unique(df_coll$experiment)){
  test = df_coll[df_coll$experiment==i & df_coll$beta > 0,]
  cr = toupper(strsplit(i, split = '_')[[1]][2])
  ls = list()
  for (j in c('conformation', 'contact', 'sandwich')){
    tf = test$gene[test$protein_change==j]
    if (length(tf)>0){
      ls[[j]] = tf
    }
  }
  if (length(ls)>1){
    print(cr)
    g = ggvenn(ls) + labs(title=cr)
    plt.list[[cr]] = g
  }
}
plt.list %>% marrangeGrob(ncol=2, nrow=2, top = '') %>%
  ggsave(file.path(exp_plt_out, 'Gene_overview_pos.pdf'),
         plot=., width=10,height=10,units='in',device='pdf',dpi=300)
### breast cancer has the largest core set

## Breast cancer =============================
breast_bg = rownames(mtx)[mtx[,'BRCA']==1]
tissue_qtl = df_coll[grep('brca', df_coll$experiment),]
co_core_contact = get_intersection_eqtl(tissue_qtl, group_col = 'protein_change')
wide = spread(co_core_contact[c('gene','protein_change','beta')], key = protein_change, value=beta)
### filter for consistent genes
cst_gene = wide$gene[which(wide$conformation * wide$contact > 0 & wide$conformation * wide$sandwich > 0)]
#### all genes are consistent
wide = subset(wide, wide$gene %in% cst_gene)
wide[['mean_beta']] = rowMeans(wide[,2:4])
spc_gene_pool = as.data.frame(table(tissue_qtl$gene))
spc_gene_pool = spc_gene_pool$Var1[which(spc_gene_pool$Freq==1)]

ont = 'BP'
#### Enrichment of the union ====
bg_pos = unique(tissue_qtl$gene[tissue_qtl$beta>0])
bg_neg = unique(tissue_qtl$gene[tissue_qtl$beta<0])
pos_GO_union = do_GO(bg_pos, background = breast_bg, ont=ont)
neg_GO_union = do_GO(bg_neg, background = breast_bg, ont=ont)
pos_GO_union_pcs = pcs_GO_out(pos_GO_union, filename='BRCA/GO_pos_union.pdf', dir=plot_out)
neg_GO_union_pcs = pcs_GO_out(neg_GO_union, filename='BRCA/GO_neg_union.pdf', dir=plot_out)

#### Core set of genes ====
pos_GO_core = do_GO(wide$gene[wide$mean_beta>0], background = bg_pos, ont=ont)
dotplot(pos_GO_core)
pos_GO_core_pcs = pcs_GO_out(pos_GO_core, filename='BRCA/GO_pos_coreSet.pdf', dir=plot_out)
neg_GO_core = do_GO(wide$gene[wide$mean_beta<0], background = bg_neg, ont=ont)
dotplot(neg_GO_core)
neg_GO_core_pcs = pcs_GO_out(neg_GO_core, filename='BRCA/GO_neg_coreSet.pdf', dir=plot_out) # nothing!

gsea_test = do_GSEA(wide, rank_nm = 'mean_beta')
gseaplot2(gsea_test, 2)
rs = gsea_test@result
p53_gs = rs[nrow(rs),'core_enrichment']
wide[wide$gene %in% strsplit(p53_gs, split='/')[[1]],]

ccle.core.pos = ccle.core$gene[which(ccle.core$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif > 0)]
ccle.contact.pos = ccle.contact$gene[which(ccle.contact$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif > 0)]
ccle.core.neg = ccle.core$gene[which(ccle.core$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif < 0)]
ccle.contact.neg = ccle.contact$gene[which(ccle.contact$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif < 0)]

ccle.tcga.co.pos = intersect(wide$gene, intersect(ccle.core.pos, ccle.contact.pos))
ccle.tcga.co.neg = intersect(wide$gene, intersect(ccle.core.neg, ccle.contact.neg))
wide = wide[wide$gene %in% c(ccle.tcga.co.neg, ccle.tcga.co.pos),]

#### Contact specific ====
contact_spc = tissue_qtl[which(tissue_qtl$protein_change == 'contact' & 
                          (!tissue_qtl$gene %in% co_core_contact$gene) &
                            (tissue_qtl$gene %in% spc_gene_pool)),]
go_con_pos = do_GO(contact_spc$gene[contact_spc$beta > 0], background = bg_pos)
go_con_neg = do_GO(contact_spc$gene[contact_spc$beta < 0], background = bg_neg) # nothing

gsea_test = do_GSEA(contact_spc[contact_spc$beta!=0,], rank_nm = 'beta')

go_con_pos_pcs = pcs_GO_out(go_con_pos, filename = 'BRCA/GO_pos_contact.pdf', dir=plot_out)
go_con_neg_pcs = pcs_GO_out(go_con_neg, filename = 'BRCA/GO_neg_contact.pdf', dir=plot_out) # nothing
dotplot(go_con_pos) # cell migration
gene_ls = unique(unlist(sapply(pcs_GO_out(go_con_pos)[[2]]$geneID, 
              function(x){return(strsplit(x, split='/')[[1]])})))
dotplot(go_con_neg) # 3 ER regulating genes, one term
cnetplot(go_con_pos)

#### core specific ====
core_spc = tissue_qtl[which(tissue_qtl$protein_change == 'conformation' & 
                            (!tissue_qtl$gene %in% co_core_contact$gene) &
                              (tissue_qtl$gene %in% spc_gene_pool)),]
go_core_pos = do_GO(core_spc$gene[core_spc$beta > 0], background = bg_pos) # nothing
go_core_neg = do_GO(core_spc$gene[core_spc$beta < 0], background = bg_neg) # nothing

gsea_test = do_GSEA(core_spc[core_spc$beta != 0,])

go_core_pos_pcs = pcs_GO_out(go_core_pos, filename = 'BRCA/GO_pos_conform.pdf', dir=plot_out)
go_core_neg_pcs = pcs_GO_out(go_core_neg, filename = 'BRCA/GO_neg_conform.pdf', dir=plot_out)
dotplot(go_core_pos) # nothing
dotplot(go_core_neg) # nothing

#### Sandwich specific ====
sw_spc = tissue_qtl[which(tissue_qtl$protein_change == 'sandwich' & 
                              (!tissue_qtl$gene %in% co_core_contact$gene) &
                            (tissue_qtl$gene %in% spc_gene_pool)),]
go_sw_pos = do_GO(sw_spc$gene[sw_spc$beta > 0], background = bg_pos)
go_sw_neg = do_GO(sw_spc$gene[sw_spc$beta < 0], background = bg_neg)

gsea_test = do_GSEA(sw_spc[sw_spc$beta != 0,])
gsea_test@result$ID
gseaplot2(gsea_test, 2)

go_sw_pos_pcs = pcs_GO_out(go_sw_pos, filename = 'BRCA/GO_pos_sanswich.pdf', dir=plot_out)
go_sw_neg_pcs = pcs_GO_out(go_sw_neg, filename = 'BRCA/GO_neg_sandwich.pdf', dir=plot_out)
dotplot(go_sw_pos)
dotplot(go_sw_neg)


save(go_sw_pos_pcs, go_sw_neg_pcs, go_core_pos_pcs, go_core_neg_pcs, 
     go_con_pos_pcs, go_con_neg_pcs, pos_GO_core_pcs, neg_GO_core_pcs,
     pos_GO_union_pcs, neg_GO_union_pcs,
     file=file.path(data_out, 'Breast_ContactCore_noCCLE_filter.RData'))

### Further analysis in breast cancer =================
load(file.path(data_out, 'Breast_ContactCore_noCCLE_filter.RData'))
# co_Gene_motile = ext_gene_GO(go_con_pos_pcs$result$geneID[1:5], do_intersect = T)
co_Gene_motile = contact_spc$gene[contact_spc$beta > 0]

brca.gene = read.table(file.path(eqtl_out, 'tcga_brca_raw_seq', 'trans_eqtl_fdr005.txt'), header = T)
brca.gene = subset(brca.gene, brca.gene$beta > 0 & brca.gene$protein_change %in% 
                     c('contact', 'conformation', 'sandwich'))

# pick genes in 1 to 6 GO terms as migration-specific signature for contact mutants
gene_in_GO = ext_gene_GO(go_con_pos_pcs$result$geneID[1:6], do_intersect = F)
brca.gene = brca.gene[brca.gene$gene %in% gene_in_GO,]
co_Gene_motile = brca.gene$gene

brca.gene = brca.gene[brca.gene$gene %in% co_Gene_motile,]
# gsea = do_GSEA(brca.gene[brca.gene$gene %in% gene_in_GO,], rank_nm = 'beta')

tcga.brca = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData', 
                            ann_bin_mut_list = c('contact', 'conformation', 'sandwich'), 
                            mode='tcga')
muts = c('has_conformation', 'has_contact', 'has_sandwich')

exp_SLC = t(assay(tcga.brca[[1]][co_Gene_motile,,'RNA']))
exp_SLC = cbind(exp_SLC, tcga.brca[[1]]@colData[rownames(exp_SLC),c('p53_state', muts)])
exp_SLC = as.data.frame(exp_SLC)
exp_SLC = subset(exp_SLC, !exp_SLC$p53_state %in% c('others', 'frameshift'))

exp_SLC$gp = exp_SLC$p53_state
for (i in muts){
  exp_SLC[which(exp_SLC[,i]==1), 'gp'] = i
}
table(exp_SLC$gp)
exp_SLC['sample'] = rownames(exp_SLC)

### Correlation between VAF and expression
vaf = tcga.brca[[2]]@data[,c('Tumor_Sample_Barcode', 'VAF', 'Hugo_Symbol')]
vaf = vaf[vaf$Hugo_Symbol=='TP53' & vaf$Tumor_Sample_Barcode %in% exp_SLC$sample]
which(duplicated(vaf$Tumor_Sample_Barcode)) # should be no duplicate samples 
exp_vaf = inner_join(exp_SLC, vaf, by=c('sample'='Tumor_Sample_Barcode'))
hist(exp_vaf$VAF) # normal

cor_mtx = matrix(NA, nrow=length(co_Gene_motile), ncol=length(unique(exp_vaf$gp)))
rownames(cor_mtx) = co_Gene_motile
colnames(cor_mtx) = unique(exp_vaf$gp)
cor_mtx = as.data.frame(cor_mtx)
pcor_mtx = cor_mtx
for (i in unique(exp_vaf$gp)){
  exp_vaf_sub = exp_vaf[exp_vaf$gp==i,]
  for (j in 1:length(co_Gene_motile)){
    tobj = cor.test(exp_vaf_sub[,j], exp_vaf_sub[,'VAF'])
    cor_mtx[co_Gene_motile[j],i] = tobj$estimate
    pcor_mtx[co_Gene_motile[j],i] = tobj$p.value
  }
}
pheatmap(cor_mtx, border_color = F, show_rownames = F)
pcor_mtx_mask = pcor_mtx
pcor_mtx_mask[pcor_mtx<0.05] = 1
pcor_mtx_mask[pcor_mtx>=0.05] = 0
pheatmap(pcor_mtx_mask, border_color = F, show_rownames = F,
         cluster_rows = T, cluster_cols = F)
colMeans(cor_mtx)
cor_mtx_long = gather(cor_mtx, key='gp', value='cor', 1:ncol(cor_mtx))
model = aov(cor~gp, cor_mtx_long)
plot(model)
TukeyHSD(model)
ggplot(gather(cor_mtx, key='gp', value='cor', 1:ncol(cor_mtx))) +
  geom_boxplot(aes(x=gp, y=cor))
test = do_GO(intersect(co_Gene_motile, rownames(cor_mtx)[order(cor_mtx$has_contact, decreasing = T)[1:100]]), 
             background = rownames(cor_mtx))
dotplot(test)

library(pheatmap)
# draw heatmap per sample
exp_SLC = exp_SLC[order(exp_SLC$gp),]
ann_row = data.frame(exp_SLC$gp, row.names = rownames(exp_SLC))
g = pheatmap(exp_SLC[,1:length(co_Gene_motile)], cluster_rows = T, cluster_cols = T, show_rownames = F,
         annotation_row = ann_row)
g = ggplotify::as.ggplot(g)
ggsave(file.path(plot_out, 'BRCA/contact_6goTerm_genes.pdf'),
       plot=g, width=8,height=6,units='in',device='pdf',dpi=300)

# draw boxplot
mycomparisons = list()
from = 'has_contact'
to = c('has_sandwich', 'has_conformation', 'Wildtype', 'nonsense')
for (i in 1:length(to)){
  mycomparisons[[i]] = c(from, to[i])
}

exp_SLC_plt = gather(exp_SLC, key='gene', value='value', 1:length(co_Gene_motile))

boxplot
g = ggplot(exp_SLC_plt, aes(x=gp, y=value)) +
  geom_boxplot() +
  facet_wrap(~gene, ncol=3) +
  stat_compare_means(comparisons = mycomparisons,
                     method = 't.test') +
  theme(axis.text.x = element_text(angle=45)) +
  mytme
g %>% ggsave(file.path(plot_out, paste('BRCA-contact-migrationGenes.pdf', sep='')),
             plot=., width=12,height=30,units='in',device='pdf',dpi=300)

exp_SLC_hm = exp_SLC_plt %>% group_by(gene, gp) %>% summarise(mean_epr = mean(value))
exp_SLC_hm = spread(exp_SLC_hm, key='gp', value='mean_epr')
exp_SLC_hm = as.data.frame(exp_SLC_hm)
rownames(exp_SLC_hm) = exp_SLC_hm$gene
exp_SLC_hm = exp_SLC_hm[,-1]

# draw heatmap mean
tqsub = subset(tissue_qtl, tissue_qtl$gene %in% rownames(exp_SLC_hm))
exp_SLC_hm = exp_SLC_hm[tqsub$gene[order(tqsub$beta, decreasing = T)],]
g = pheatmap(exp_SLC_hm, labels_row = rownames(exp_SLC_hm), 
             #labels_col = '', 
             cluster_cols = F, cluster_rows = F, border_color = F,
             width = 6, height = 10)
g = ggplotify::as.ggplot(g)
# g = g + scale_x_discrete(labels = colnames(exp_SLC_hm))
# filename = file.path(plot_out, 'contact_GO1_2_genes.pdf'
ggsave(file.path(plot_out, 'BRCA/contact_GO1_2_genes.pdf'),
       plot=g, width=4,height=4,units='in',device='pdf',dpi=300)
write.table(exp_SLC_hm %>% mutate(gene=rownames(exp_SLC_hm)), 
            file.path(data_out, 'BRCA_contact_spc_GOI.txt'), sep='\t', row.names = F, quote = F)


