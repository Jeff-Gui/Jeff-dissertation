library(tidyverse)
library(ggpubr)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('enrich_utils.R')
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'GO')
data_out = file.path(dir_home, 'data_out')

flt_df = function(go_df, cutoff=0.01){
  bg = sapply(go_df$GeneRatio, function(x){eval(parse(text = x))})
  go_df = go_df[which(bg > cutoff),]
  return(go_df)
}

## eQTL mapping outputs

### Experiment: contrast samples with p53 mutation to those without.
# does not filter beta here
beta_cutoff = 0.5
coll = load_eQTL_output(eqtl_out, beta=beta_cutoff, exclude = 'tcga_nine_pool')
study = c()
map_mut = c()
num_mut = c()
total_gene = c()
df_coll = data.frame()
for (i in names(coll)){
  nm = toupper(strsplit(i, split='_')[[1]][2])
  study = c(study, nm)
  df = coll[[i]]
  if (class(df)=='logical'){
    map_mut = c(map_mut, NA)
    num_mut = c(num_mut, NA)
    total_gene = c(total_gene, NA)
  } else {
    df = subset(df, abs(df$beta) > beta_cutoff)
    df_coll = rbind(df_coll, df)
    map_mut = c(map_mut, paste(unique(df$protein_change), collapse = '|'))
    num_mut = c(num_mut, length(unique(df$protein_change)))
    total_gene = c(total_gene, nrow(df))
  }
}
sry = data.frame('Study'=study, 'Mapped_mutant'=map_mut, 'Num_mutant'=num_mut, 
                 'Num_gene' = total_gene)
sry[['Avg_gene']] = sry$Num_gene / sry$Num_mutant
print(sry)


sry = sry[which(!is.na(sry$Study)),]
sry = sry[which(sry$Study != 'NINE'),]
df_coll = df_coll[grep('tcga', df_coll$experiment),]
if (length(grep('nine', df_coll$experiment))>0){
  df_coll = df_coll[-grep('nine', df_coll$experiment),]
}
write.table(sry, file.path(dir_home, 'output_summary.txt'), row.names = F, quote=F, sep='\t')


## do volcano plot
hist(df_coll$beta,breaks=100)
ggplot(df_coll) +
  geom_histogram(aes(x=beta), binwidth = 0.01)
g = ggplot(df_coll) +
  geom_point(aes(x=beta, y=-log10(FDR), color=experiment),size=0.5) +
  geom_vline(xintercept = 0.5, linetype='dotted', color='black') +
  mytme + theme(legend.text = element_text(size=12)) +
  geom_vline(xintercept = -0.5, linetype='dotted', color='black')
ggsave(file.path(dirname(plot_out), 'Volcano_overview.png'),
       plot=g, width=8.27*1.2,height=11.69*0.5,units='in',device='png',dpi=300)

df_coll_hs = subset(df_coll, df_coll$protein_change=='hot_spot')
ggplot(df_coll) +
  # geom_hline(yintercept = -log10(0.05), linetype='dotted', color='grey') +
  geom_point(aes(x=beta, y=-log10(FDR), color=protein_change)) +
  geom_vline(xintercept = 0.5, linetype='dotted', color='black') +
  geom_vline(xintercept = -0.5, linetype='dotted', color='black') +
  facet_wrap(~experiment, scale='free') +
  mytme + theme(
    legend.text = element_text(size=10),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=14, color='black', face='bold'),
    panel.border = element_rect(color='black', size=1, fill='transparent')
    #legend.position = 'none'
  )

## Gene ontology enrichment and also enrichment of transcription factors

### Run GO
cache = FALSE
if (cache){
  load(file.path(dir_home, 'GO_result.RData'))
  # load(file.path(dir_home, 'GO_result_no_BG_filter.RData'))
} else {
  term2gene = load_TRUST_term2gene()
  beta_cutoff = beta_cutoff # analyse genes with beta value above the cutoff
  min_gene = 20  # TODO?
  result = data.frame()
  go_coll = list()
  for (epr in unique(df_coll$experiment)){
    df = subset(df_coll, df_coll$experiment == epr)
    for (mut in unique(df$protein_change)){
      print(paste(epr, mut, sep='_'))
      # num of pos/neg correlated genes and enriched GO terms
      msg = rep(0,6)
      df_sub = subset(df, df$protein_change == mut)
      df_sub = subset(df_sub, abs(df_sub$beta)> beta_cutoff)
      df_up = subset(df_sub, beta > beta_cutoff)
      df_down = subset(df_sub, beta < -beta_cutoff)
      msg[1] = nrow(df_up)
      msg[2] = nrow(df_down)
      cot = 3
      fsx = c('_pos_GO.pdf', '_pos_TF_GO.pdf', '_neg_GO.pdf', '_neg_TF_GO.pdf')
      fsx_shot = c('pos', 'pos_TF', 'neg', 'neg_TF')
      for (gp in list(df_up, df_down)){
        if (nrow(gp)>min_gene){
          up_GO = do_GO(gp)
          po = pcs_GO_out(up_GO, pvaluecutoff = 0.05, 
                          filename = paste(epr, '_', mut, fsx[cot-2], sep=''), 
                          dir = plot_out)
          msg[cot] = po$n_sig_term
          if (msg[cot]>0){
            go_coll[[paste(epr, mut, fsx_shot[cot-2], sep='-')]] = po$result
          }
        }
        if (nrow(gp)>min_gene){
          if (length(intersect(gp$gene, term2gene$gene))!=0){
            ern = enricher(gp$gene, TERM2GENE = term2gene, pvalueCutoff = 0.05)
            if (is.null(ern)){
              msg[cot+1] = 0
            } else {
              po = pcs_GO_out(ern, pvaluecutoff = 0.05, 
                              filename = paste(epr, '_', mut, fsx[cot-1], sep=''), 
                              dir = plot_out)
              msg[cot+1] = po$n_sig_term
            }
          } else {
            msg[cot+1] = 0
          }
          if (msg[cot+1]>0){
            go_coll[[paste(epr, mut, fsx_shot[cot-1], sep='-')]] = po$result
          }
        }
        cot = cot + 2
      }
      result = rbind(result, c(epr, mut, msg))
    }
  }
  colnames(result) = c('Experiment', 'Mutation', 'Pos_gene', 'Neg_gene', 'Pos_GO', 'Pos_TF_GO', 'Neg_GO', 'Neg_TF_GO')
  gc()
  result[,3:8] = as.numeric(as.matrix(result[,3:8]))
  save(go_coll, result, file = file.path(dir_home, 'GO_result_no_BG_filter.RData'))
  
  ## QC gene ontology, require 5% gene ratio
  to_remove = c()
  sign_map = c('pos','neg','pos_TF','neg_TF')
  names(sign_map) = c('Pos_GO', 'Neg_GO', 'Pos_TF_GO', 'Neg_TF_GO')
  for (i in 1:nrow(result)){
    for (j in c('Pos_GO', 'Neg_GO', 'Pos_TF_GO', 'Neg_TF_GO')){
      if (result[[j]][i] > 0){
        sign = sign_map[j]
        cd = paste(result$Experiment[i], result$Mutation[i], sign, sep='-')
        go_coll[[cd]] = flt_df(go_coll[[cd]], cutoff = 0.05)
        result[[j]][i] = nrow(go_coll[[cd]])
        if (nrow(go_coll[[cd]])==0){
          to_remove = c(to_remove, cd)
        }
      }
    }
  }
  if (!is.null(to_remove)){
    go_coll = go_coll[which(!names(go_coll) %in% to_remove)]
  }
  save(go_coll, result, file = file.path(dir_home, 'GO_result.RData'))
}


load(file.path(dir_home, 'GO_result_no_BG_filter.RData'))
# plot eQTL result
# print(result)
result['code'] = sapply(result$Experiment, function(x){
  return(toupper(strsplit(x, split = '_')[[1]][2]))
})
result_back = result
result = result_back
exclude = c()
# exclude = c('HNSC', 'STAD', 'LUAD') # too few mutations
name_map = read.csv('/Users/jefft/Desktop/p53_project/datasets/mut_name_map.csv', header=T)
mutation_to_plt = name_map$label
names(mutation_to_plt) = name_map$code
result = subset(result, !result$code %in% exclude)
result = subset(result, result$Mutation %in% names(mutation_to_plt))
result$Mutation = mutation_to_plt[result$Mutation]
result_long = gather(result, key='gp', value='count', 3:8)
result_long = result_long[-grep('TF', result_long$gp),]
result_long$count = as.numeric(result_long$count)
result_long$gp = factor(result_long$gp, levels = c('Pos_gene', 'Neg_gene', 'Pos_GO', 'Neg_GO'))
to_plt = c('Pos_gene', 'Neg_gene')
to_plt = c('Pos_GO', 'Neg_GO')
strc_mut_plt = c('DNA contact', 'Core')
other_mut_plt = mutation_to_plt[which(!mutation_to_plt %in% c(cons_mut_plt, strc_mut_plt))]

plt.list = list()
count = 1

for (mut_plt in list(strc_mut_plt, other_mut_plt)){
  for (to_plt in list(c('Pos_gene', 'Neg_gene'), c('Pos_GO', 'Neg_GO'))){
    if (length(grep('gene', to_plt))>0){
      lb = 'Gene'
    } else {
      lb = 'GO'
    }
    g = ggplot(subset(result_long, gp %in% to_plt & Mutation %in% mut_plt),
               aes(x=Mutation, y=count, fill=gp)) +
      geom_bar(stat = 'identity', position = 'dodge') +
      facet_wrap(~code, scales = 'free') +
      scale_fill_manual(values = c('coral', 'royalblue'), name = '', 
                        labels = c('Positive', 'Negative')) +
      labs(x='', y='Count', title = lb) +
      mytme +
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            strip.background = element_blank(),
            strip.text = element_text(face='bold', size=14),
            panel.border = element_rect(color='black', size=1, fill=NA)) +
      scale_y_log10(expand = c(0,0))
    plt.list[[count]] = g
    count = count + 1
  }
}

plt.list %>% marrangeGrob(ncol=1, nrow=2, top = '') %>%
  ggsave(file.path(dirname(plot_out), 'Overview.pdf'),
         plot=., width=11.69,height=8.27*2,units='in',device='pdf',dpi=300)


## Overlapping of hot spot mutants
# get overlap among cancers
# set 0.5 beta cutoff
df_coll = subset(df_coll, abs(df_coll$beta) > 0.5)
library(ggvenn)
exps = c('tcga_blca_raw_seq', 
         'tcga_brca_raw_seq', 'tcga_lgg_raw_seq', 'tcga_coad_raw_seq')
hs = subset(df_coll, df_coll$protein_change=='hot_spot' & 
              df_coll$experiment %in% exps)
coll_pos_hs = list()
for (i in exps){
  nm = toupper(strsplit(i, split = '_')[[1]][2])
  coll_pos_hs[[nm]] = hs$gene[which(hs$beta > 0 & hs$experiment == i)]
}
ggvenn(coll_pos_hs, stroke_color = 'white') %>%
  ggsave(file.path(dirname(plot_out), 'Gene_overlap_hotspot_pos.pdf'),
         plot=., width=6,height=5,units='in',device='pdf',dpi=300)



## print overlap among cancers
# exclude stomach cancer that only has 31 genes
hs_inter = get_intersection_eqtl(hs, group_col = 'experiment', 
                                 doPlot = F, check_beta = 'pos')
ups = paste(unique(hs_inter$gene[which(hs_inter$beta > 0)]), collapse = ', ')
print(paste('Positive related genes in BLCA, BRCA, COAD, LGG by hot spot mutations:', ups))
hs_inter = get_intersection_eqtl(hs, group_col = 'experiment', 
                                 doPlot = F, check_beta = 'neg')
downs = paste(unique(hs_inter$gene[which(hs_inter$beta < 0)]), collapse = ', ')
print(paste('Negative related genes in BLCA, BRCA, COAD, LGG by hot spot mutations:', downs))


# Positively correlated genes:
#   
# - CLCN2: "Voltage-gated chloride channel". __Sig in all 4, but not same sign.__
# 
# - ESYT3: "Sphingolipid metabolism (REACTOME) and Metabolism. Gene Ontology (GO) annotations related to this gene include lipid binding and phospholipid binding".
# 
# - HDGF: "A member of the hepatoma-derived growth factor family. The encoded protein has mitogenic and DNA-binding activity and may play a role in cellular proliferation and differentiation. High levels of expression of this gene enhance the growth of many tumors". __Sig in all 4, same positive sign, but not significant in CCLE.__
# 
# - NOL10: "Polyglutamine Binding Protein 5".
# 
# - UBXN2A: "Ubiquitin binding and acetylcholine receptor".

### Test if they are operational in CCLE.

# load CCLE
### Check in CCLE, load data
load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
# Annotate p53 state
if (!'p53_state' %in% colnames(dt[[1]]@colData)){
  p53_ann = annotate_sample_mut(dt[[2]]@data)
  dt[[1]]@colData[['p53_state']] = 'Wildtype'
  for (i in names(p53_ann)){
    dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
  }
}
meta_mut_home = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_mut_fp = list.files(meta_mut_home)
meta_mut_fp = sapply(meta_mut_fp, 
                     function(x){return(file.path(meta_mut_home, x))})
meta_mut = load_meta_mut(meta_mut_fp)



# Get genes to test
breast_tcga_hot_spot_pos = subset(df_coll, experiment=='tcga_brca_raw_seq' &
                                    abs(beta) > 0.5 & protein_change == 'hot_spot')
genes = breast_tcga_hot_spot_pos$gene




# run t-test in CCLE
check_site = 'Breast'
check_mut = 'hot_spot'
# genes = c('EDA2R', 'MDM2', 'SPATA18')

idx = rownames(dt[[1]]@colData)
if (!is.null(check_site)){
  idx_ftd = idx[which(dt[[1]]@colData$PRIMARY_SITE %in% check_site)]
} else {
  idx_ftd = idx
}
if (length(grep('p\\.', check_mut))==0){
  mut_aa = unique(meta_mut$aa_pos[meta_mut$meta_mut_id == check_mut])
  b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                  snp_list = list(mut_aa),
                                  samples = idx_ftd, 
                                  mode = 'position')
} else {
  b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                  snp_list = list(check_mut),
                                  snp_col_nm = 'Protein_Change',
                                  samples = idx_ftd, 
                                  mode = 'amino_acid')
}
mutant = rownames(b_m)[which(b_m==1)]
wildtype = intersect(idx_ftd, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'Wildtype')])
nonsense = intersect(idx_ftd, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'nonsense')])

# check which specific mutation mut+ cells have
# mut_state = sapply(names(b_m_ftd)[which(b_m_ftd=='mut+')], function(x){
#   sub_tb = subset(dt[[2]]@data, dt[[2]]@data$Tumor_Sample_Barcode==x &
#            dt[[2]]@data$Hugo_Symbol=='TP53')
#   return(sub_tb$Protein_Change)
# })

count = 0
gene_pool = t(assay(dt[[1]][genes,c(wildtype, mutant, nonsense),'RNA']))
ps = c()
fcs = c()
ps_null = c()
fcs_null = c()
for (i in 1:length(genes)){
  if (!genes[i] %in% colnames(gene_pool)){
    ps = c(ps, NA)
    fcs = c(fcs, NA)
    ps_null = c(ps_null, NA)
    fcs_null = c(fcs_null, NA)
    count = count + 1
    print(paste(genes[i], 'not in expression matrix.'))
  } else {
    a = gene_pool[wildtype, genes[i]]
    b = gene_pool[mutant, genes[i]]
    c = gene_pool[nonsense, genes[i]]
    p = t.test(a,b)$p.value
    ps = c(ps, p)
    fc = mean(b) - mean(a)
    fcs = c(fcs, fc)
    ps_null = c(ps_null, t.test(b,c)$p.value)
    fcs_null = c(fcs_null, mean(b) - mean(c))
  }
}
ccle.valid = data.frame('gene' = genes, 'ccle.wt.p' = ps, 'ccle.wt.dif' = fcs, 
                        'ccle.wt.FDR' = p.adjust(ps, method = 'fdr'),
                        'ccle.null.p' = ps_null, 'ccle.null.dif' = fcs_null, 
                        'ccle.null.FDR' = p.adjust(ps_null, method = 'fdr'))
print(count)
ccle.valid_signf = subset(ccle.valid, ccle.wt.FDR < 0.05)
ccle.valid_signf_null = subset(ccle.valid, ccle.wt.FDR < 0.05 & ccle.null.p < 0.05)



### integrate CCLE and TCGA result
rownames(ccle.valid_signf) = ccle.valid_signf$gene
rownames(breast_tcga_hot_spot_pos) = breast_tcga_hot_spot_pos$gene
itg = cbind(ccle.valid_signf[,-which(colnames(ccle.valid_signf)=='gene')],
            breast_tcga_hot_spot_pos[rownames(ccle.valid_signf),])
itg = subset(itg, itg$ccle.wt.dif * itg$beta > 0 &
               itg$ccle.null.dif * itg$ccle.wt.dif > 0)
itg[['GO_ann']] = ann_GO(itg$gene)

itg_pos = subset(itg, beta>0)
test = do_GO(itg_pos)
g = dotplot(test)
ggsave(file.path(dirname(plot_out), 'TCGA_BRCA-hot_spot-pos-CCLE_psig.pdf'),
       width=6, height=5, units='in', device='pdf', dpi=300, plot = g)
rs = test@result
rs = subset(rs, rs$p.adjust < 0.05)
enriched_genes = unique(unlist(sapply(rs$geneID, function(x){return(strsplit(x, split='/')[[1]])})))
gc()
write.table(itg_pos %>% mutate(GO_ann = as.character(GO_ann)),
            file.path(data_out, 'TCGA_BRCA-hot_spot-pos-CCLE_sig.txt'),
            sep='\t', quote=F, row.names=F)
itg_pos = read.table(file.path(data_out, 'TCGA_BRCA-hot_spot-pos-CCLE_sig.txt'), sep='\t', header = T)

itg_neg = subset(itg, beta<0)
test = do_GO(itg_neg)
g = dotplot(test)
# ggsave(file.path(dirname(plot_out), 'TCGA_BRCA-hot_spot-neg-CCLE_null_sig.pdf'),
#        width=6, height=5, units='in', device='pdf', dpi=300, plot = g)
# rs = test@result
# rs = subset(rs, rs$p.adjust < 0.05)
# enriched_genes = unique(unlist(sapply(rs$geneID, function(x){return(strsplit(x, split='/')[[1]])})))
# gc()
write.table(itg_neg %>% mutate(GO_ann = as.character(GO_ann)),
            file.path(data_out, 'TCGA_BRCA-hot_spot-neg-CCLE_sig.txt'),
            sep='\t', quote=F, row.names=F)

## Test RNAi
rnai = load_rnai()
cell_pool = list(wildtype, mutant, nonsense)
for (i in 1:3){
  cell_pool[[i]] = intersect(cell_pool[[i]], colnames(rnai))
}
names(cell_pool) = c('wildtype', 'mutant', 'nonsense')

gene_to_check = intersect(itg_pos$gene,rownames(rnai)) # 1/133 checked missing
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

### GO annotation, RNAi, deprecated
# sub_idt_sig[['GO_ann']] = ann_GO(sub_idt_sig$gene)

# # annotate null mutation
# p53_null_ccle = get_null_samples(dt[[2]]@data)
# 
# 
# plt.list = list()
# gene_done = c()
# for (i in 1:nrow(sub_idt_sig)){
#   if (sub_idt_sig$gene[i] %in% gene_done){next}
#   gene_done = c(gene_done, sub_idt_sig$gene[i])
#   plt_df = cbind(t(assay(dt[[1]][sub_idt_sig$gene[i],idx_ftd,'RNA'])), 
#                  b_m[idx_ftd,])
#   plt_df = as.data.frame(plt_df)
#   if (ncol(plt_df)==1){
#     print(paste(sub_idt_sig$gene[i], 'missing from the expression set.'))
#     next
#   }
#   colnames(plt_df) = c('Gene', 'Mut')
#   plt_df$Gene = as.numeric(plt_df$Gene)
#   plt_df_gn = cbind(dt[[1]]@colData[idx_ftd,], plt_df)
#   plt_df_gn = as.data.frame(plt_df_gn)
#   plt_df_gn$cell_nm = sapply(rownames(plt_df_gn), function(x)return(strsplit(x, split='_')[[1]][1]))
#   # plt_df_gn[['null']] = FALSE
#   plt_df_gn[p53_null_ccle, 'Mut'] = 'null'
#   g = ggplot(plt_df_gn, aes(x=Mut, y=Gene)) + theme_classic() +
#     geom_violin() +
#     stat_compare_means(comparisons = list(c('mut-','mut+'),
#                                           c('null', 'mut+')), 
#                        method = 't.test') +
#     geom_jitter(size=0.5, width=0.2) +
#     geom_boxplot(width=0.2, outlier.shape = NA) +
#     labs(x=check_mut, y=sub_idt_sig$gene[i]) +
#     #geom_text(aes(label=cell_nm, alpha=plt_df_gn$Mut=='mut+'), size=2, nudge_x=0, nudge_y=0.02, angle=30) +
#     scale_alpha_manual(values = c(0,1)) + 
#     theme(legend.position = 'none')
#   plt.list[[sub_idt_sig$gene[i]]] = g
# }
# 
# marrangeGrob(grobs=plt.list,ncol=3,nrow=3) %>% 
#   ggsave(file.path(plot_out, 'TCGA_Pan_hot_spot_CCLE_Pan_test.pdf'),
#          plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)

## Test the trend in CCLE RNAi data

# gene_tar = 'CLCN2'
# cl_n = paste(mut_aa, collapse = '_')
# rnai_plt = as.data.frame(dt[[1]]@colData[idx_ftd,])
# rnai_plt['mutation_state'] = b_m[idx_ftd, cl_n]
# rnai_plt[p53_null_ccle,'mutation_state'] = 'null'
# 
# rnai = load_rnai()
# co_cell = intersect(idx_ftd, colnames(rnai))
# rnai_plt = rnai_plt[co_cell,]
# rnai_plt = cbind(rnai_plt, t(rnai[gene_tar, rownames(rnai_plt)]))
# colnames(rnai_plt)[which(colnames(rnai_plt) == gene_tar)] = 'Gene'
# rnai_plt = rnai_plt[-which(is.na(rnai_plt$Gene)),]
# print(table(rnai_plt$mutation_state))
# 
# g = ggplot(rnai_plt, aes(x=mutation_state, y=Gene)) + theme_classic() +
#   geom_violin() +
#   stat_compare_means(comparisons = list(c('mut-','mut+'),
#                                         c('null', 'mut+')), 
#                      method = 't.test') +
#   geom_jitter(size=0.5, width=0.2) +
#   geom_boxplot(width=0.2, outlier.shape = NA) +
#   labs(x=check_mut, y=gene_tar, title = 'CCLE RNAi') +
#   scale_alpha_manual(values = c(0,1)) + 
#   # facet_wrap(~PRIMARY_SITE) +
#   theme(legend.position = 'none')
# g
# ggsave(file.path(plot_out, 'CLCN2_CCLE_rnai_pan.pdf'),
#        plot=g, width=6,height=5,units='in',device='pdf',dpi=300)
# inspect = subset(rnai_plt, rnai_plt$mutation_state=='mut+')
# inspect = inspect[order(inspect$Gene, decreasing = T),c('Gene', 'PRIMARY_SITE')]
