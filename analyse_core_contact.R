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
library(ggsci)
library(ComplexUpset)
source('../overlap_utils.R')
exp_plt_out = file.path(dir_home, 'plots', 'coreVScontact')

# GO ====
gl = load_GO_out(exp_home = dir_home, gene_ratio_min = 0,
                 fname = 'GO_BP_result_no_BG_filter.RData')
go_df = gl$df
go_df = go_df[grep('GO', go_df$ID),]
signn = 'neg'
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

## GO enrichment for core mutation found in most cancers ====
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

# Gene ====
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
sign = 'neg'
plt.list = list()
hs_mtx_coll = list()
cancer_OI = c('BRCA','COAD', 'LGG', 'STAD')
for (i in unique(df_coll$experiment)){
  if (sign == 'pos'){
    test = df_coll[df_coll$experiment==i & df_coll$beta > 0,]
  } else {
    test = df_coll[df_coll$experiment==i & df_coll$beta < 0,]
  }
  cr = toupper(strsplit(i, split = '_')[[1]][2])
  if (!cr %in% cancer_OI){next}
  ls = list()
  for (j in c('conformation', 'contact', 'sandwich')){
    tf = test$gene[test$protein_change==j]
    if (length(tf)>0){
      ls[[j]] = tf
    }
  }
  hs_mtx_coll[[cr]] = gen_upSet_mtx(ls)
  if (length(ls)>1){
    print(cr)
    g = ggvenn(ls, stroke_color = 'white') + labs(title=cr) + scale_fill_d3(palette = "category20")
    plt.list[[cr]] = g
  }
}
save(hs_mtx_coll, file = file.path(data_out, paste('UpsetMtx_strc_',sign,'.RData', sep='')))
plt.list %>% marrangeGrob(ncol=2, nrow=2, top = '') %>%
  ggsave(file.path(exp_plt_out, paste('Gene_overview_', sign, '.pdf', sep='')),
         plot=., width=10,height=10,units='in',device='pdf',dpi=300)

## Find significant overlaps ====
ovl_p_coll = list()
for (i in 1:length(hs_mtx_coll)){
  if (ncol(hs_mtx_coll[[i]])==3){
    ovl_p = run_overlap_test(hs_mtx_coll[[i]], n_loop = 1000, test_condition = c(T,T,T))
    ovl_p_coll[[names(hs_mtx_coll)[i]]] = 1-ovl_p$p.left
  }
}
for (i in 1:length(ovl_p_coll)){
  if (ovl_p_coll[[i]]<0.05){
    print(names(ovl_p_coll)[[i]])
  }
}
# Positive gene overlap: BRCA, COAD, LGG, STAD
# Negative gene overlap: BRCA, COAD, LGG, STAD

# [[checkpoint]] ====
load(file.path(data_out, 'UpsetMtx_strc_neg.RData'))
neg_um = hs_mtx_coll
load(file.path(data_out, 'UpsetMtx_strc_pos.RData'))
pos_um = hs_mtx_coll
remove(hs_mtx_coll)

# Breast cancer ====
ont = 'BP'
sub_folder = 'LGG' # BRCA

uni_bg = rownames(mtx)[mtx[, sub_folder]==1]
tissue_qtl = df_coll[grep(tolower(sub_folder), df_coll$experiment),]
co_core_contact = get_intersection_eqtl(tissue_qtl, group_col = 'protein_change')
wide = spread(co_core_contact[c('gene','protein_change','beta')], key = protein_change, value=beta)

#### Filter consistent genes ====
cst_gene = wide$gene[which(wide$conformation * wide$contact > 0 & wide$conformation * wide$sandwich > 0)]
print(paste('All in the core set:', nrow(wide), ', Same sign:', length(cst_gene)))
uncon_wide = wide[!wide$gene %in% cst_gene,]
print(uncon_wide) # NFATC1, PDZD4
wide = subset(wide, wide$gene %in% cst_gene)
wide[['mean_beta']] = rowMeans(wide[,2:4])

#### Enrichment of the union ====
pos_um_tis = pos_um[[sub_folder]]
neg_um_tis = neg_um[[sub_folder]]
bg_pos = rownames(pos_um_tis)
bg_neg = rownames(neg_um_tis)
pos_GO_union = do_GO(bg_pos, background = uni_bg, ont=ont)
neg_GO_union = do_GO(bg_neg, background = uni_bg, ont=ont)
pos_GO_union_pcs = pcs_GO_out(pos_GO_union, filename=file.path(sub_folder, 'GO_pos_union.pdf'), dir=plot_out, cneplt=F)
neg_GO_union_pcs = pcs_GO_out(neg_GO_union, filename=file.path(sub_folder, 'GO_neg_union.pdf'), dir=plot_out, cneplt=F)

#### Specific set of genes ====
ern_df = data.frame()
colnames(pos_um_tis)
test_cdts = list(c(T,T,T), c(T,F,F), c(F,T,F), c(F,F,T))
file_prefix = c('Core', 'Conform', 'Contact', 'Sandwich')
ont = 'BP'
for (i in 1:length(test_cdts)){
  pos_spc = bg_pos[get_idx_condt(pos_um_tis,test_condition = test_cdts[[i]])]
  neg_spc = bg_neg[get_idx_condt(neg_um_tis,test_condition = test_cdts[[i]])]
  pos_GO = do_GO(pos_spc, background = bg_pos, ont=ont)
  neg_GO = do_GO(neg_spc, background = bg_neg, ont=ont)
  pos_GO_pcs = pcs_GO_out(pos_GO, filename=file.path(sub_folder, paste('GO_', ont, '_pos_', file_prefix[i], '.pdf', sep='')), dir=plot_out, cneplt = F)
  neg_GO_pcs = pcs_GO_out(neg_GO, filename=file.path(sub_folder, paste('GO_', ont, '_neg_', file_prefix[i], '.pdf', sep='')), dir=plot_out, cneplt = F)
  ern_df = rbind(ern_df, c(pos_GO_pcs$n_sig_term, neg_GO_pcs$n_sig_term, length(pos_spc), length(neg_spc)))
}
colnames(ern_df) = c('Pos_GO', 'Neg_GO', 'Pos_gene', 'Neg_gene')
ern_df$experiment = file_prefix
ern_df = rbind(ern_df, c(pos_GO_union_pcs$n_sig_term, neg_GO_union_pcs$n_sig_term, length(bg_pos), length(bg_neg), 'Union'))

write.table(ern_df,file.path(plot_out, sub_folder, 'summary_enrichment.txt'), sep='\t', quote = F, row.names = F)

# Further analysis in breast cancer ==========
## Plot upset ====
upset(pos_um_tis, intersect = colnames(pos_um_tis),
      min_degree = 1, mode='exclusive_intersection')

pos_um_tis_ctr = load_controls(um_to_merge = pos_um_tis, ng_ctrs = T)
### where are the controls
pgs = pos_um_tis_ctr[,c("conformation","contact","sandwich", "pg")]
upset(pgs, intersect = colnames(pgs),
      mode = 'exclusive_intersection')
neg_um_tis_ctr = load_controls(um_to_merge = neg_um_tis, ng_ctrs = T)
### where are the controls
ngs = neg_um_tis_ctr[,c("conformation","contact","sandwich", "ng1", "ng2")]
upset(ngs, intersect = colnames(ngs),
      mode = 'exclusive_intersection')
colSums(ngs)


colnames(pos_um_tis_ctr)
pos_um_tis_ctr = pos_um_tis_ctr[get_idx_condt(pos_um_tis_ctr, c(NA,NA,NA,F,F,NA,NA,NA,NA,F,T,F)),]
test = pos_um_tis_ctr[pos_um_tis_ctr$pg==1,]

colnames(pos_um_tis)
test = do_GO(rownames(pos_um_tis)[get_idx_condt(pos_um_tis, c(F,T,F))], background = bg_pos)
dotplot(test)
## pick adhesion and migration related terms
gene_in_GO = ext_gene_GO(test@result$geneID[c(3,4,6,7)], do_intersect = F)
pos_um_tis_mig = pos_um_tis_ctr[rownames(pos_um_tis_ctr) %in% gene_in_GO,]
colnames(pos_um_tis_mig)
mig_gene_eqtl = df_coll[df_coll$gene %in% gene_in_GO & df_coll$experiment=='tcga_brca_raw_seq',]
mig_gene_eqtl = cbind(mig_gene_eqtl[,1:7], pos_um_tis_mig[mig_gene_eqtl$gene,])
write.table(mig_gene_eqtl, file.path(data_out, 'BRCA_migration_qtl.txt'), sep='\t', quote = F, row.names = F)
# overlap with R273H peaks
GOI = rownames(pos_um_tis_mig)[pos_um_tis_mig$peak.over.chek1==1]

brca.gene = read.table(file.path(eqtl_out, 'tcga_brca_raw_seq', 'trans_eqtl_fdr005.txt'), header = T)
brca.gene = subset(brca.gene, brca.gene$beta > 0 & brca.gene$protein_change %in% 
                     c('contact', 'conformation', 'sandwich'))
brca.gene[brca.gene$gene %in% GOI,]

colnames(pos_um_tis)
contact_pos_genes = bg_pos[get_idx_condt(pos_um_tis,
                          test_condition = c(F,T,F))]
contact_pos_GO = do_GO(contact_pos_genes, background = bg_pos, ont='BP')
sandwich_neg_genes = bg_neg[get_idx_condt(neg_um_tis,
                            test_condition = c(F,T,F))]
sandwich_neg_GO = do_GO(sandwich_neg_genes, background = bg_neg, ont='BP')

## Profile migraiton gene signature ====
mut_groups = c('contact', 'conformation', 'sandwich')
dt.tcga = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData', 
                          ann_bin_mut_list = mut_groups, mode='tcga')
dt.ccle = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd/BRCA/clean_data.RData',
                          ann_bin_mut_list = mut_groups, mode = 'ccle')
df.tcga = cbind(t(assay(dt.tcga[[1]][GOI,,'RNA'])), as.data.frame(dt.tcga[[1]]@colData[,c('p53_state', 'has_contact', 'has_conformation', 'has_sandwich')]))
df.ccle = cbind(t(assay(dt.ccle[[1]][GOI,,'RNA'])), as.data.frame(dt.ccle[[1]]@colData[,c('p53_state', 'has_contact', 'has_conformation', 'has_sandwich')]))
# !!! choose one !!!
df.plt = df.ccle
df.plt = df.tcga
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
# CCLE width: 5, TCGA width: 10
ggsave(file.path(plot_out,'BRCA', 'HeatmapTCGA.pdf'),
       plot = g + theme(legend.key.width = unit(0.4,units = 'inch')),
       width=8.27*2, height=6, units='in', device='pdf', dpi=300, bg = 'transparent')
ggsave(file.path(plot_out,'BRCA', 'WntHeatmapCCLE.pdf'),
       plot = g + theme(legend.key.width = unit(0.3,units = 'inch')),
       width=8.27*0.6, height=4, units='in', device='pdf', dpi=300, bg = 'transparent')

