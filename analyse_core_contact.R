setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'coreVScontact')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)

gl = load_GO_out(exp_home = dir_home, gene_ratio_min = 0)
go_df = gl$df
go_df = go_df[grep('GO', go_df$ID),]

#### Core VS Contact mutant, GO =====================================
signn = 'neg'
exp_plt_out = file.path(dir_home, 'plots', 'coreVScontact')

cc_go_df = subset(go_df, go_df$mutation %in% c('conformation', 'contact', 'sandwich'))
cc_go_df = cc_go_df[grep('GO', cc_go_df$ID),]

cc_go_df = subset(cc_go_df, cc_go_df$sign==signn)
table(cc_go_df$experiment)

plt.list = list()
for (i in unique(cc_go_df$cancer)){
  test = cc_go_df[cc_go_df$cancer==i,]
  g = ggvenn(list('core'=test$ID[test$mutation=='core'], 
                  'contact'=test$ID[test$mutation=='contact'])) +
    labs(title=i)
  plt.list[[i]] = g
}
plt.list %>% marrangeGrob(ncol=3, nrow=3, top = '') %>%
  ggsave(file.path(exp_plt_out, 'GO_overview_neg.pdf'),
         plot=., width=10,height=10,units='in',device='pdf',dpi=300)

# GO enrichment for core mutation found in most cancers
signn = 'neg'
core_gos = go_df[go_df$mutation == 'core',]
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
g %>% ggsave(file.path(exp_plt_out, paste(signn, '_pan_core_GO_dotPlot.pdf', sep='')),
             plot=., width=25,height=15,units='in',device='pdf',dpi=300)


# There is only GO enrichment of contact mutant in BRCA and LGG.
# load CCLE
load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
ccle.core = test_ccle(dt, check_mut = 'core', genes=NULL, check_site = 'Breast')
ccle.contact = test_ccle(dt, check_mut = 'contact', genes=NULL, check_site = 'Breast')
ccle.core = ccle.core$ccle.valid
ccle.contact = ccle.contact$ccle.valid
save(ccle.core, ccle.contact, file = file.path(data_out, 'ccle_breast_ContactCorevsWT_allgenes.RData'))

### Gene ================
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

## Gene overlap overview
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

## Genes in breast cancer
breast_bg = rownames(mtx)[mtx[,'BRCA']==1]
tissue_qtl = df_coll[grep('brca', df_coll$experiment),]
co_core_contact = get_intersection_eqtl(tissue_qtl, group_col = 'protein_change')
wide = spread(co_core_contact[c('gene','protein_change','beta')], key = protein_change, value=beta)
### filter for consistent genes
cst_gene = wide$gene[which(wide$conformation * wide$contact > 0 & wide$conformation * wide$sandwich > 0)]
#### all genes are consistent
wide = subset(wide, wide$gene %in% cst_gene)
wide[['mean_beta']] = rowMeans(wide[,2:4])

#------------ enrichment of the union
bg_pos = unique(tissue_qtl$gene[tissue_qtl$beta>0])
bg_neg = unique(tissue_qtl$gene[tissue_qtl$beta<0])
pos_GO_union = do_GO(bg_pos, background = breast_bg)
neg_GO_union = do_GO(bg_neg, background = breast_bg)
pcs_GO_out(pos_GO_union, filename='GO_pos_union.pdf', dir=plot_out)
pcs_GO_out(neg_GO_union, filename='GO_neg_union.pdf', dir=plot_out)

#------------ the overlapping core set of genes
pos_GO_core = do_GO(wide$gene[wide$mean_beta>0], background = bg_pos)
dotplot(pos_GO_core)
pcs_GO_out(pos_GO_core, filename='GO_pos_coreSet.pdf', dir=plot_out)
neg_GO_core = do_GO(wide$gene[wide$mean_beta<0], background = bg_neg)
dotplot(neg_GO_core)
pcs_GO_out(neg_GO_core, filename='GO_neg_coreSet.pdf', dir=plot_out)


ccle.core.pos = ccle.core$gene[which(ccle.core$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif > 0)]
ccle.contact.pos = ccle.contact$gene[which(ccle.contact$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif > 0)]
ccle.core.neg = ccle.core$gene[which(ccle.core$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif < 0)]
ccle.contact.neg = ccle.contact$gene[which(ccle.contact$ccle.wt.p<0.05 & ccle.core$ccle.wt.dif < 0)]

ccle.tcga.co.pos = intersect(wide$gene, intersect(ccle.core.pos, ccle.contact.pos))
ccle.tcga.co.neg = intersect(wide$gene, intersect(ccle.core.neg, ccle.contact.neg))
wide = wide[wide$gene %in% c(ccle.tcga.co.neg, ccle.tcga.co.pos),]

#------------ contact specific
contact_spc = tissue_qtl[which(tissue_qtl$protein_change == 'contact' & 
                          (!tissue_qtl$gene %in% co_core_contact$gene)),]
go_con_pos = do_GO(contact_spc$gene[contact_spc$beta > 0], background = bg_pos)
go_con_neg = do_GO(contact_spc$gene[contact_spc$beta < 0], background = bg_neg) # nothing

# go_con_pos = do_GO(intersect(contact_spc$gene[contact_spc$beta > 0], 
#                       ccle.contact.pos))
# go_con_neg = do_GO(intersect(contact_spc$gene[contact_spc$beta < 0], 
#                              ccle.contact.neg))

pcs_GO_out(go_con_pos, filename = 'GO_pos_contact.pdf', dir=plot_out)
pcs_GO_out(go_con_neg, filename = 'GO_neg_contact.pdf', dir=plot_out)
dotplot(go_con_pos) # cell migration
gene_ls = unique(unlist(sapply(pcs_GO_out(go_con_pos)[[2]]$geneID, 
              function(x){return(strsplit(x, split='/')[[1]])})))
dotplot(go_con_neg) # 3 ER regulating genes, one term
cnetplot(go_con_pos)

#-------- core specific
core_spc = tissue_qtl[which(tissue_qtl$protein_change == 'conformation' & 
                            (!tissue_qtl$gene %in% co_core_contact$gene)),]
go_core_pos = do_GO(core_spc$gene[core_spc$beta > 0], background = bg_pos) # nothing
go_core_neg = do_GO(core_spc$gene[core_spc$beta < 0], background = bg_neg) # nothing

# go_core_pos = do_GO(intersect(core_spc$gene[core_spc$beta > 0], 
#                              ccle.core.pos))
# go_core_neg = do_GO(intersect(core_spc$gene[core_spc$beta < 0], 
#                              ccle.core.neg))

pcs_GO_out(go_core_pos, filename = 'GO_pos_conform.pdf', dir=plot_out)
pcs_GO_out(go_core_neg, filename = 'GO_neg_conform.pdf', dir=plot_out)
dotplot(go_core_pos) # nothing
dotplot(go_core_neg) # nothing
# cnetplot(go_core_pos)

#-------- sandwich specific
sw_spc = tissue_qtl[which(tissue_qtl$protein_change == 'sandwich' & 
                              (!tissue_qtl$gene %in% co_core_contact$gene)),]
go_sw_pos = do_GO(sw_spc$gene[sw_spc$beta > 0], background = bg_pos) # nothing
go_sw_neg = do_GO(sw_spc$gene[sw_spc$beta < 0], background = bg_neg) # nothing

pcs_GO_out(go_sw_pos, filename = 'GO_pos_sanswich.pdf', dir=plot_out)
pcs_GO_out(go_sw_neg, filename = 'GO_neg_sandwich.pdf', dir=plot_out)
dotplot(go_sw_pos) # nothing
dotplot(go_sw_neg) # nothing


save(go_con_pos, go_con_neg, go_core_pos, go_core_neg, file=file.path(data_out, 'Breast_ContactCore_noCCLE_filter.RData'))

load(file.path(data_out, 'Breast_ContactCore.RData'))
strsplit(go_con_pos@result$geneID[1], split='/')[[1]]

