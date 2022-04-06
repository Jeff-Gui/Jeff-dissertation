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
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)
library(ggrepel)

# Load Gene ========
beta_cutoff = 0
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
df_coll[['cancer']] = sapply(df_coll$experiment, function(x){
  return(strsplit(x, split='_')[[1]][2])
})
df_coll$cancer = toupper(df_coll$cancer)

# Load GO ========
load(file.path(dir_home, 'GO_BP_result_no_BG_filter.RData'))
go_df = data.frame()
for (i in 1:length(go_coll)){
  nm = names(go_coll)[i]
  sub_df = go_coll[[nm]]
  sub_df[['experiment']] = nm
  rownames(sub_df) = NULL
  go_df = rbind(go_df, sub_df)
}
hist(go_df$Count, breaks=100)
go_df['cancer'] = sapply(go_df$experiment, 
                         function(x){toupper(strsplit(x, split = '_')[[1]][2])})
go_df['sign'] = 'neg'
go_df[grep('pos', go_df$experiment), 'sign'] = 'pos'
go_df['mutation'] = sapply(go_df$experiment,
                           function(x){strsplit(x, split='-')[[1]][2]})

# Overlapping gene ========
library(UpSetR) # upset plot requires sample-group matrix
# https://blog.csdn.net/tuanzide5233/article/details/83109527
library(ggvenn)
#exps = unique(df_coll$experiment)
exps = c('BLCA', 'STAD', 'BRCA', 'LGG', 'COAD')
hs = subset(df_coll, df_coll$protein_change=='hot_spot' & 
              df_coll$cancer %in% exps)
upset_mtx = matrix(0, nrow=length(unique(hs$gene)), ncol=length(exps))
rownames(upset_mtx) = unique(hs$gene)
colnames(upset_mtx) = exps
upset_mtx = as.data.frame(upset_mtx)
coll_pos_hs = list()
for (i in exps){
  hits = hs$gene[which(hs$beta < 0 & hs$cancer == i)]
  upset_mtx[hits, i] = 1
  coll_pos_hs[[i]] = hits
}
upset_mtx = upset_mtx[-which(rowSums(upset_mtx)==0),]

ggvenn(coll_pos_hs, stroke_color = 'white') %>%
  ggsave(file.path(plot_out, 'Gene_overlap_hotspot_neg.pdf'),
         plot=., width=6,height=5,units='in',device='pdf',dpi=300)

pdf(file=file.path(plot_out, "Gene_overlap_US_hotspot_neg.pdf"), onefile=FALSE, width = 18, height = 6)
upset(upset_mtx, text.scale = c(2,2,2,2,2,2), mb.ratio = c(0.6,0.4))
dev.off()

## print overlap among cancers
hs_inter = get_intersection_eqtl(hs, group_col = 'experiment', 
                                 doPlot = F, check_beta = 'pos')
ups = paste(unique(hs_inter$gene[which(hs_inter$beta > 0)]), collapse = ', ')
print(paste('Positive related genes in BLCA, BRCA, COAD, LGG by hot spot mutations:', ups))
hs_inter = get_intersection_eqtl(hs, group_col = 'experiment', 
                                 doPlot = F, check_beta = 'neg')
downs = paste(unique(hs_inter$gene[which(hs_inter$beta < 0)]), collapse = ', ')
print(paste('Negative related genes in BLCA, BRCA, COAD, LGG by hotspot mutations:', downs))

test = do_GO(unique(hs_inter$gene[which(hs_inter$beta > 0)]), 
             background = hs$gene[which(hs$beta>0)])


### Test hot spot genes in BRCA in CCLE ========
tcga.brca = subset(df_coll, df_coll$protein_change=='hot_spot' & df_coll$experiment=='tcga_brca_raw_seq')

#### load control genes
library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
ctrs = ctrs[-grep('wt', ctrs$Gene_annotation),]
pos.ctrs = ctrs$Gene

#### load TCGA-BRCA hotspot VS nonsense (p sig, no p adj)
tcga.brca.ns = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-null/outputs/tcga_brca_raw_seq/trans_eqtl.txt', header = T)
tcga.brca.ns = tcga.brca.ns[tcga.brca.ns$protein_change=='hot_spot',]
tcga.ns.up = tcga.brca.ns$gene[tcga.brca.ns$beta > 0]
#### load pre-computed CCLE BRCA test
ccle.brca = read.table(file.path(CCLE_home, 'BRCA/hotspot_comp.txt'), sep='\t', header = T)
#### load literature control
known_target = read.table('/Users/jefft/Desktop/p53_project/datasets/Fischer2017/TableS2.txt', header = T)
#### load H1299 R273H ChIP-seq
peaks = read.table('/Users/jefft/Desktop/p53_project/datasets/PRJEB20314/p53-1-2-merged_Peaks.txt', header = T)

known_wt_up = known_target$Gene.Symbol[known_target$sum.direct.regulation.score...16..0..16. >= 2]
known_wt_down = known_target$Gene.Symbol[known_target$sum.direct.regulation.score...16..0..16. <= -2]
mean_chek1 = mean(peaks$mean_ern.rep1[peaks$gene=='CHEK1'], peaks$mean_ern.rep2[peaks$gene=='CHEK1'])
peaks = na.omit(peaks)
peak_over_chek1 = peaks$gene[rowMeans(as.matrix(peaks[,c('mean_ern.rep1', 'mean_ern.rep2')])) >= mean_chek1]
tcga.up = tcga.brca$gene[tcga.brca$beta > 0]
ccle.brca.mask = ccle.brca
ccle.brca.mask[is.na(ccle.brca.mask)] = 0
ccle.rna.up = ccle.brca.mask$gene[ccle.brca.mask$rna_hs.wt_dif > 0 & ccle.brca.mask$rna_hs.wt_p < 0.05]
ccle.rnai.down = ccle.brca.mask$gene[ccle.brca.mask$rnai_hs.wt_dif < 0 & ccle.brca.mask$rnai_hs.wt_p < 0.05]
ccle.crispr.down = ccle.brca.mask$gene[ccle.brca.mask$crispr_hs.wt_dif < 0 & ccle.brca.mask$crispr_hs.wt_p < 0.05] # nothing if use padj
ccle.rna.up.null = ccle.brca.mask$gene[ccle.brca.mask$rna_hs.ns_dif > 0 & ccle.brca.mask$rna_hs.ns_p < 0.05]

mylist = list('pos.ctrs' = pos.ctrs,
              'tcga.up' = tcga.up, 
              # 'tcga.ns.up' = tcga.ns.up,
              'ccle.rna.up' = ccle.rna.up, 'ccle.rna.ns.up' = ccle.rna.up.null,
              'ccle.rnai.down' = ccle.rnai.down, 'ccle.crispr.down' = ccle.crispr.down,
              'known_wt_up' = known_wt_up, 'known_wt_down' = known_wt_down,
              'peaks_over_chek1' = peak_over_chek1)

library(ComplexUpset)
um = gen_upSet_mtx(mylist)
um = um[um$known_wt_up==0 & um$known_wt_down==0,]
um = um[,-which(colnames(um) %in% c('known_wt_up', 'known_wt_down'))]
print(colSums(um))

#pdf(file=file.path(plot_out, 'upset_brca.pdf'), width=12, height=6)
#dev.off()

# https://krassowski.github.io/complex-upset/articles/Examples_R.html
g = upset(um[rowSums(um)>=1,], colnames(um), min_size = 1, width_ratio = 0.2,
      sort_intersections_by = 'degree', sort_intersections = 'ascending',
      queries = list(upset_query(set='pos.ctrs', color='coral', fill='coral'),
                     upset_query(set='tcga.up', color='royalblue', fill='royalblue'),
                     upset_query(intersect=c('tcga.up', 'peaks_over_chek1'), fill='orange', color='orange'),
                     upset_query(intersect=c('tcga.up', 'ccle.rna.up'), fill='orange', color='orange'),
                     upset_query(intersect=c('tcga.up', 'ccle.rna.up', 'ccle.rna.ns.up'), fill='orange', color='orange'),
                     upset_query(intersect=c('tcga.up', 'peaks_over_chek1', 'ccle.rna.up'), fill='orange', color='orange')),
      base_annotations=list(
        'Intersection size'=intersection_size(
          mode = 'inclusive_intersection',
          text=list(vjust=-0.1,hjust=0,angle=45,size=3),
          #mapping=aes(fill=mpaa),
          text_colors=c(on_background='brown', on_bar='coral'))
          + annotate(
            geom='text', x=Inf, y=Inf,
            label=paste('Total genes:', nrow(um)),
            vjust=1, hjust=1
          ))
      )
ggsave(file.path(plot_out, 'upset.pdf'),
       plot=g, width=11.69*1.2,height=8.27*0.7,units='in',device='pdf',dpi=300)


### Shuffle upset matrix, how many observed is significant? ====
print(colnames(um))
# T: must be 1, F: must be 0, NA: no matter. 
test_condition = c(F, T, T, NA, NA, NA, NA)
overlap_stat = run_overlap_test(um, test_condition = test_condition, n_loop = 1000)

### Enrichment analysis ====
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
bg = rownames(mtx)[which(mtx$BRCA==1)]
# TCGA RNA (vs WT), CCLE RNA (vs WT), H1299 ChIP, CCLE RNA (vs Null)
# qua_gene = rownames(um)[get_idx_condt(um, test_condition = c(T,T,T,NA,NA,T))] # 21 genes
# TCGA RNA (vs WT), CCLE RNA (vs WT), H1299 ChIP
trip_gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,T,T,NA,NA,NA,T))] # 96 genes
trip_gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,T,T,T,NA,NA,NA))] # 100 genes
# TCGA RNA and H1299 ChIP (no CCLE)
dual_gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,T,NA,NA,NA,NA,T))] # 620 genes
dual_gene = rownames(um)[get_idx_condt(um, test_condition = c(NA,T,T,NA,NA,NA,NA))] # 547 genes

# test = do_GO(qua_gene, background = bg, ont='BP')
test = do_GO(trip_gene, background = bg, ont='BP')
dotplot(test)
test@result$geneID[1]
test = do_GO(dual_gene, background = bg)
g = dotplot(test)
ggsave(file.path(plot_out, 'TCGA_BRCA-hot_spot-pos-CCLE_psig_GO.pdf'),
       width=8, height=8, units='in', device='pdf', dpi=300, plot = g)

tcga.brca_test = subset(tcga.brca, tcga.brca$gene %in% trip_gene)





# Archived ====
# load CCLE
load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
p53_ann = annotate_sample_mut(dt[[2]]@data)
dt[[1]]@colData[['p53_state']] = 'Wildtype'
for (i in names(p53_ann)){
  dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
}

# generate ccle valid hotspot VS wildtype
check_site = 'Breast'
check_mut = 'hot_spot'
ccle.valid.all = test_ccle(dt, check_mut, genes=NULL, check_site=check_site)
write.table(ccle.valid.all[[1]], file.path(data_out, 'ccle_breast_hotspotVSwt_TCGA-BRCA-genes.txt'), sep='\t', row.names = F, quote = F)
ccle.valid.all = ccle.valid.all[[1]]
# mutant: 10, wildtype: 19, nonsense: 4 for breast cancer cell lines

# Get genes to test
breast_tcga_hot_spot_pos = subset(df_coll, experiment=='tcga_brca_raw_seq' &
                                    abs(beta) > 0 & protein_change == 'hot_spot')
genes = breast_tcga_hot_spot_pos$gene
# genes = c('EDA2R', 'MDM2', 'SPATA18')

# run t-test in CCLE
check_site = 'Breast'
check_mut = 'hot_spot'

ccle.valid.brca = test_ccle(dt, check_mut, genes=genes, check_site=check_site)
print(length(ccle.valid.brca$unmatched))
ccle.valid.brca = ccle.valid.brca$ccle.valid
#### !!! change file name according to experiment
write.table(ccle.valid, file.path(data_out, 'ccle_breast_hotspotVSwt_TCGA-BRCA-genes.txt'), sep='\t', row.names=F, quote=F)

# ccle.valid.tcga = read.table(file.path(data_out, 'ccle_breast_hotspotVSwt_TCGA-BRCA-genes.txt'), header = T)

## Load CCLE valid data =============================
ccle.valid.all = read.table(file.path(ccle_home, 'ccle_breast_hotspotVSwt_allgenes.txt'), header = T)
ccle.valid.all[['found.in.tcga']] = F
ccle.valid.all$found.in.tcga[which(ccle.valid.all$gene %in% tcga.brca$gene)] = T

# we need to filter those found in both, but with different sign...
co_genes = ccle.valid.all$gene[ccle.valid.all$found.in.tcga]
idx_ist = which(deframe(ccle.valid.all[which(ccle.valid.all$gene %in% co_genes), c('gene','ccle.wt.dif')])[co_genes] * 
  deframe(tcga.brca[which(tcga.brca$gene %in% co_genes),c('gene', 'beta')])[co_genes] < 0)
print(length(idx_ist))
ccle.valid.all['wrong_sign'] = FALSE
ccle.valid.all[which(ccle.valid.all$gene %in% co_genes[idx_ist]), 'wrong_sign'] = T
# 1923 of them 
ccle.valid.all = ccle.valid.all[!ccle.valid.all$wrong_sign,]

table(ccle.valid.all$found.in.tcga)
ggplot(ccle.valid.all) +
  geom_point(aes(x=ccle.wt.dif, y=-log10(ccle.wt.FDR), color=found.in.tcga), size=0.5) +
  geom_hline(yintercept = -log10(0.05), color='red')
cross_tb = table(ccle.valid.all$ccle.wt.FDR < 0.05, ccle.valid.all$found.in.tcga)
rownames(cross_tb) = c('NOT_FDR_SIG', 'FDR_SIG')
colnames(cross_tb) = c('NOT_IN_TCGA', 'IN_TCGA')
chisq.test(cross_tb)
print(cross_tb)
# Considering wrong sign mut VS wt
#              NOT_IN_TCGA  IN_TCGA (59 not detected in CCLE)
# NOT_FDR_SIG        11575     5789
# FDR_SIG               88      188
# considering wrong sign
#              NOT_IN_TCGA  IN_TCGA (59 not detected in CCLE)
# NOT_FDR_SIG        14363     3803
# FDR_SIG              128      148

odd = (cross_tb[2,2] / cross_tb[1,2]) / (cross_tb[2,1] / cross_tb[1,1])
print(odd) # 4.27

# ccle.valid_signf = subset(ccle.valid, ccle.wt.FDR < 0.05)
# ccle.valid_signf_null = subset(ccle.valid, ccle.wt.FDR < 0.05 & ccle.null.p < 0.05)

co_signif = subset(ccle.valid.all, ccle.valid.all$ccle.wt.FDR<0.05 & ccle.valid.all$found.in.tcga & 
                     abs(ccle.valid.all$ccle.wt.dif) > 0 & !ccle.valid.all$wrong_sign)
# test = do_GSEA(co_signif, rank_nm = 'ccle.wt.dif')

### integrate CCLE and TCGA result
bg = unique(c(ccle.valid.all$gene[which(ccle.valid.all$ccle.wt.FDR<0.05)], tcga.brca$gene))
itg = inner_join(co_signif, tcga.brca, by=c('gene'='gene'))
# itg = subset(itg, itg$ccle.wt.dif * itg$beta > 0 &
#               itg$ccle.null.dif * itg$ccle.wt.dif > 0)
# itg[['GO_ann']] = ann_GO(itg$gene)

itg_pos = subset(itg, beta>0)
test = do_GO(itg_pos, background = bg)
g = dotplot(test)
ggsave(file.path(dirname(plot_out), 'TCGA_BRCA-hot_spot-pos-CCLE_psig.pdf'),
       width=6, height=5, units='in', device='pdf', dpi=300, plot = g)
write.table(itg_pos, file.path(data_out, 'TCGA_BRCA-hot_spot-pos-CCLE_sig.txt'),
            sep='\t', quote=F, row.names=F)

itg_neg = subset(itg, beta<0)
test = do_GO(itg_neg, background = bg)
g = dotplot(test) # nothing
ggsave(file.path(dirname(plot_out),'hotspot', 'TCGA_BRCA-hot_spot-neg-CCLE_psig.pdf'),
       width=6, height=5, units='in', device='pdf', dpi=300, plot = g)
write.table(itg, file.path(data_out, 'TCGA_BRCA-hot_spot_CCLE_sig.txt'), sep='\t', quote=F, row.names = F)


## genes significant when comparing to null
itg_sig_null = subset(itg, itg$ccle.null.p < 0.05 & (itg$ccle.wt.dif * itg$ccle.null.dif > 0))
test = do_GO(itg_sig_null$gene[itg_sig_null$beta > 0], background=itg$gene[itg$beta>0])
itg_sig_null$sig_mul = itg_sig_null$ccle.null.p * itg_sig_null$FDR
itg_sig_null$dif_mul = abs(itg_sig_null$ccle.null.dif) * abs(itg_sig_null$beta) * ifelse(itg_sig_null$beta > 0, 1,-1)
itg_sig_null_to_plt = itg_sig_null %>% mutate(lb = gene)
itg_sig_null_to_plt$lb[abs(itg_sig_null$dif_mul)<=1] = NA
write.table(itg_sig_null_to_plt, 
            file.path(data_out, 'TCGA_BRCA-hot_spot_CCLE_null_psig.txt'), sep='\t', quote=F, row.names = F)
g = ggplot(itg_sig_null_to_plt,aes(x=dif_mul, y=-log10(sig_mul))) +
  geom_point() +
  geom_text_repel(aes(label=lb)) +
  labs(x = 'Beta (eQTL) * Difference (DE)', y= '-log10(eQTL_FDR * DE_p)', 
       title =str_wrap('Significant genes Mut:WT in TCGA (eQTL) and Mut:Nonsense in CCLE (DE)', width = 50)) +
  mytme

ggsave(file.path(dirname(plot_out), 'hotspot', 'TCGA_BRCA-hot_spot-pos-CCLE_null_psig.pdf'),
       width=8, height=5, units='in', device='pdf', dpi=300, plot = g)



# Test RNAi
breast_cells = rownames(dt[[1]]@colData)[which(dt[[1]]@colData$PRIMARY_SITE=='Breast')]
rnai = load_rnai()
itg_pos = read.table(file.path(data_out, 'TCGA_BRCA-hot_spot-pos-CCLE_sig.txt'), sep='\t', header = T)
b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, snp_list = list(c(175,245,248,249,273,282)),
                                samples = breast_cells, 
                                mode = 'position')
mutant = rownames(b_m)[which(b_m==1)]

wildtype = intersect(breast_cells, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'Wildtype')])
nonsense = intersect(breast_cells, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'nonsense')])
cell_pool = list(wildtype, mutant, nonsense)
for (i in 1:3){
  cell_pool[[i]] = intersect(cell_pool[[i]], colnames(rnai))
}
names(cell_pool) = c('wildtype', 'mutant', 'nonsense')

gene_to_check = intersect(itg_pos$gene,rownames(rnai)) # 6/73 checked missing
gene_pool_rnai = rnai[gene_to_check,unlist(cell_pool)]
for (i in 1:nrow(gene_pool_rnai)){
  a = as.numeric(na.omit(gene_pool_rnai[i,cell_pool$wildtype]))
  b = as.numeric(na.omit(gene_pool_rnai[i,cell_pool$mutant]))
  if (!is.na(mean(a)) & !is.na(mean(b))){
    if (mean(a) > mean(b)){
      p = t.test(a,b)$p.value
      if (p < 0.05){
        print(rownames(gene_pool_rnai)[i])
        print(p)
      }
    }
  }
}
# [1] "KLF11"
#  0.0384017
#  "PTPRK"
#  0.02341799
#  "CEBPB"
#  0.001786938

## Overlap GO among cancers ====================

### Intersect among cancers
gl = load_GO_out(exp_home = exp_home, gene_ratio_min = 0)
go_df = gl$df

plot_out = file.path(exp_home, 'plots', 'hotspot')
# take intersect
library(ggvenn)
signn = 'pos'
hot_spot_go = go_df[grep('hot_spot', go_df$experiment),]
hot_spot_go = subset(hot_spot_go, hot_spot_go$sign==signn)
hot_spot_go = hot_spot_go[grep('GO', hot_spot_go$ID),]
coll = list()
for (i in unique(hot_spot_go$cancer)){
  coll[[i]] = hot_spot_go$ID[which(hot_spot_go$cancer == i)]
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, paste(signn, '_pan_GO_overlap.pdf', sep='')),
             plot=., width=5,height=5,units='in',device='pdf',dpi=300)

#### dot plot of each top 5 terms for each group
signn = 'pos'
top_go_plt = data.frame()
for (i in c('BLCA', 'COAD', 'LGG', 'BRCA')){
  for (j in c('hot_spot')){ # 'core', 'contact'
    sub_go_df = subset(go_df, (go_df$cancer == i) & (go_df$mutation == j) & (go_df$sign == signn))
    sub_go_df = sub_go_df[grep('GO', sub_go_df$ID),]
    if (nrow(sub_go_df)>0){
      top_go_plt = rbind(top_go_plt, sub_go_df[order(sub_go_df$p.adjust)[1:(min(5,nrow(sub_go_df)))],])
    }
  }
}

# top_go_plt$Description = factor(top_go_plt$Description, levels = unique(top_go_plt$Description))
g = ggplot(top_go_plt, aes(x=mutation)) +
  geom_point(aes(y=Description, color=-log10(p.adjust), size=gene_ratio_val, group=cancer)) +
  facet_wrap(~cancer, scale='free', nrow = 2) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=25)) +
  theme(axis.text.x = element_text(angle=0, vjust = 1, size=16),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size=16, face='bold')) +
  mytme
g %>% ggsave(file.path(plot_out, paste(signn, '_pan_hotspot_GO_dotPlot.pdf', sep='')),
             plot=., width=10,height=12,units='in',device='pdf',dpi=300)

ov = c('BRCA', 'COAD')
ov = c('BRCA', 'BLCA')
coll_brca_coad = coll[which(names(coll) %in% ov)]
coll_brca_coad = Reduce(intersect, coll_brca_coad) # nothing
test = subset(hot_spot_go, ID %in% coll_brca_coad & cancer %in% ov)
test %>% group_by(ID) %>% summarize(mn = mean(gene_ratio_val), sd=sd(gene_ratio_val))
g = ggplot(test, aes(x=cancer)) +
  geom_point(aes(y=Description, color=-log10(p.adjust), size=gene_ratio_val, group=cancer)) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=35)) +
  theme(axis.text.x = element_text(angle=0, vjust = 1, size=16),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size=16, face='bold')) +
  mytme
signn = 'pos'
g %>% ggsave(file.path(plot_out, paste(paste(ov, collapse = '-'), 
                                       signn, '_overlap_GO_dotPlot.pdf', sep='')),
             plot=., width=8,height=10,units='in',device='pdf',dpi=300)

## A heatmap of pvalue GO terms
mtx_pvalue = matrix(NA, nrow=length(unique(hot_spot_go$ID)), ncol=length(unique(hot_spot_go$cancer)))
rownames(mtx_pvalue) = unique(hot_spot_go$ID)
colnames(mtx_pvalue) = unique(hot_spot_go$cancer)
mtx_gnRatio = mtx_pvalue
for (i in 1:nrow(hot_spot_go)){
  mtx_pvalue[hot_spot_go$ID[i], hot_spot_go$cancer[i]] = hot_spot_go$p.adjust[i]
  mtx_gnRatio[hot_spot_go$ID[i], hot_spot_go$cancer[i]] = hot_spot_go$gene_ratio_val[i]
}
mtx_pvalue[is.na(mtx_pvalue)] = 1
mtx_gnRatio[is.na(mtx_gnRatio)] = 0
mtx_itg = -log10(mtx_pvalue) * mtx_gnRatio
pheatmap::pheatmap(mtx_itg)
pheatmap::pheatmap(mtx_pvalue)
