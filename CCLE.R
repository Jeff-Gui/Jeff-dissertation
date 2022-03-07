setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')
source('enrich_utils.R')
library(tidyverse)
library(ggpubr)
library(gridExtra)
config_name = 'ccle_inspect_only.yaml'
refresh_log = FALSE
save_intermediate = FALSE  # may spend extra time
use_cache = FALSE
use_cache_geno_pca = TRUE
source = TRUE

## Handle config
default_cfg = yaml.load_file(file.path('config', 'default.yaml'))
if (config_name != 'dafault.yaml'){
  config = yaml.load_file(file.path('config', config_name))
  config = merge_cfg(default_cfg, config)
}
config = pcs_cfg(config)
dt_cfg = config$dataset
eqtl_cfg = config$eQTL
preprocess_cfg = config$preprocess

# Load data
if (!use_cache){
  dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                     exp_nm = dt_cfg$exp_nm, 
                     cna_nm = dt_cfg$cna_nm, 
                     mut_nm = dt_cfg$mut_nm,
                     sample_meta_nm = dt_cfg$sample_meta_nm,
                     patient_meta_nm = dt_cfg$patient_meta_nm,
                     case_complete_nm = dt_cfg$case_complete_nm,
                     case_list_dir_nm = dt_cfg$case_list_dir_nm,
                     na.str = dt_cfg$na.str,
                     gene_col_nm = dt_cfg$gene_col_nm,
                     normalize_genes = preprocess_cfg$norm_gene,
                     quantile_norm = preprocess_cfg$quantile,
                     z_score = preprocess_cfg$z_score,
                     has_log_ed = dt_cfg$has_log_ed,
                     rm_low_expr_gene = as.numeric(preprocess_cfg$rm_low_expr_gene),
                     # diag_out = config$output,
                     filter_protein_coding = eqtl_cfg$genepos,
                     diploid_norm = preprocess_cfg$diploid_norm)
  gc()
  if (save_intermediate){
    save(dt, file = file.path(dirname(dt_cfg$dataset_home), 'clean_data_inspect.RData'))
  }
} else {
  load(file = file.path(dirname(dt_cfg$dataset_home), 'clean_data_inspect.RData'))
}

p53_ann = annotate_sample_mut(dt[[2]]@data)
dt[[1]]@colData[['p53_state']] = 'wildtype'
for (i in names(p53_ann)){
  dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
}

## Check positive controls
### Load controls
library(xlsx)
ctrs = read.xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx', sheetIndex = 1)
ctrs = na.omit(ctrs)
gene_pool = unique(ctrs$Gene)
mis = 100*(1 - sum(gene_pool %in% rownames(assay(dt[[1]][,,'RNA']))) / length(gene_pool))
print(paste(as.character(mis), '% missing in the expression matrix.', sep=''))
gene_pool = gene_pool[which(gene_pool %in% rownames(assay(dt[[1]][,,'RNA'])))]
gene_pool = c(gene_pool, c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1'))
ctr_expr = assay(dt[[1]][,,'RNA'])[c(unique(gene_pool), 'GPI'),]
# Combine expression with colData
ctr_expr = as.data.frame(cbind(t(ctr_expr), dt[[1]]@colData))

b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                c(as.list(paste('p.', unique(ctrs$Mutant), sep='')),
                                  list(c('p.R273H', 'p.R248Q'),
                                       c('p.R175H', 'p.Y163C', 'p.R249S', 'p.K164E', 'p.E286K', 'p.G245S'))), 
                                'Protein_Change',
                                rownames(ctr_expr))
for (i in 1:ncol(b_m)){
  b_m[,i] = as.factor(b_m[,i])
}
ctr_expr = cbind(ctr_expr, b_m)

mut = 'p.R175H_p.Y163C_p.R249S_p.K164E_p.E286K_p.G245S'
mut = 'p.R273H_p.R248Q'
gn = 'TEAD4'

ctrs = data.frame()
for (mut in c('p.R175H_p.Y163C_p.R249S_p.K164E_p.E286K_p.G245S',
              'p.R273H_p.R248Q')){
  for (gn in c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1')){
    ctrs = rbind(ctrs, c(mut, gn)) 
  }
}
colnames(ctrs) = c('Mutant', 'Gene')



plt.list = list()
tissue = sapply(rownames(ctr_expr), function(x){
  a = strsplit(x, split='_')
  return(a[[1]][length(a[[1]])])
})
cl_name = sapply(rownames(ctr_expr), function(x){
  a = strsplit(x, split='_')
  return(a[[1]][1])
})
ctr_expr['tissue'] = tissue
ctr_expr['cl_name'] = cl_name

ctr_expr = subset(ctr_expr, ctr_expr)
ctr_expr = subset(ctr_expr, ctr_expr$cl_name %in% c('MDA-MB-468', 'U373MG', 'U251MG', 'SF295', 'HCC193', 'PC-9',
                                                    'HCC1395', 'HCC1954', 'BT-549', 'LN229', 'M059J', 'M059K', 'U138MG', 'SK-MEL-2', 'SK-LMS-1') | 
                    ctr_expr$p53_state=='wildtype')

for (i in 1:nrow(ctrs)){
  # mut = paste('p.', ctrs$Mutant[i], sep='')
  mut = ctrs$Mutant[i]
  gn = ctrs$Gene[i]
  if (!gn %in% colnames(ctr_expr) | !mut %in% colnames(ctr_expr)){
    next
  }
  sub_ctr_expr = subset(ctr_expr, ctr_expr[[mut]]==2 | ctr_expr[['p53_state']] == 'wildtype')
  plt_dt = sub_ctr_expr[,c(mut, gn)]
  colnames(plt_dt) = c('Mut', 'Gene')
  plt_dt$Mut[which(plt_dt$Mut==1)] = 'wt'
  plt_dt$Mut[which(plt_dt$Mut==2)] = 'mut'
  g = ggplot(plt_dt, aes(x=Mut, y=Gene)) + theme_classic() +
    geom_violin() +
    stat_compare_means(comparisons = list(c('wt','mut'))) +
    geom_jitter(size=0.5, width=0.2) +
    geom_boxplot(width=0.2, outlier.shape = NA) +
    labs(x=mut, y=gn)
  # g
  if (!is.null(g)){
    plt.list[[as.character(i)]] = g
  }
}
marrangeGrob(grobs=plt.list,ncol=3,nrow=3) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/eQTL', 'CCLE_YAP_controls-paper.pdf'),
         plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)

  

stat_compare_means(comparisons = list(c('missense', 'wildtype'),
                                        c('missense', 'nonsense'),
                                        c('nonsense', 'wildtype')))

## Check identified eQTLs.
out_home = '/Users/jefft/Desktop/p53_project/scripts/eQTL/outputs/'
meta_mut_home = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_mut_fp = list.files(meta_mut_home)
meta_mut_fp = sapply(meta_mut_fp, function(x){return(file.path(meta_mut_home, x))})
meta_mut = load_meta_mut(meta_mut_fp)

idt_qtl = read.table(file.path(out_home, 'tcga_brca_raw_seq/trans_eqtl_fdr005.txt'),
                     header = T)
idx = rownames(dt[[1]]@colData)
site = sapply(idx, function(x)strsplit(x, split='_')[[1]][2])

check_mut = 'hot_spot'
idx_ftd = idx[site=='BREAST']
if (length(grep('p\\.', check_mut))==0){
  mut_aa = meta_mut$aa_pos[meta_mut$meta_mut_id == check_mut]
  b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                  snp_list = list(mut_aa),
                                  samples = rownames(dt[[1]]@colData), 
                                  mode = 'position')
} else {
  b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                  snp_list = list(check_mut),
                                  snp_col_nm = 'Protein_Change',
                                  samples = rownames(dt[[1]]@colData), 
                                  mode = 'amino_acid')
}
b_m[which(b_m==0)] = 'mut-'
b_m[which(b_m==1)] = 'mut+'

sub_idt = subset(idt_qtl, idt_qtl$protein_change==check_mut)
sub_idt = subset(sub_idt, abs(sub_idt$beta)>1) # filter for effect size over 1
sub_idt = sub_idt[order(sub_idt$beta, decreasing = T),]

gene_pool = t(assay(dt[[1]][unique(sub_idt$gene),idx_ftd,'RNA']))
b_m_ftd = b_m[idx_ftd,]

# check which specifci mutation mut+ cells have
sapply(names(b_m_ftd)[which(b_m_ftd=='mut+')], function(x){
  sub_tb = subset(dt[[2]]@data, dt[[2]]@data$Tumor_Sample_Barcode==x &
           dt[[2]]@data$Hugo_Symbol=='TP53')
  return(sub_tb$Protein_Change)
})

ps = c()
fcs = c()
for (i in 1:nrow(sub_idt)){
  if (!sub_idt$gene[i] %in% colnames(gene_pool)){
    ps = c(ps, NA)
    fcs = c(fcs, NA)
    print(paste(sub_idt$gene[i], 'not in expression matrix.'))
  } else {
    a = gene_pool[b_m_ftd=='mut-',sub_idt$gene[i]]
    b = gene_pool[b_m_ftd=='mut+',sub_idt$gene[i]]
    p = t.test(a,b)$p.value
    ps = c(ps, p)
    fc = mean(b) - mean(a)
    fcs = c(fcs, fc)
  }
}
sub_idt['ccle.p'] = ps
sub_idt['ccle.dif'] = fcs
sub_idt['ccle.FDR'] = p.adjust(sub_idt$ccle.p, method = 'fdr')
# dual-sig: FDR all <0.05; differene same sign
sub_idt_sig = subset(sub_idt, (ccle.FDR < 0.05) & (sub_idt$beta * sub_idt$ccle.dif > 0))
sub_idt_sig[['GO_ann']] = ann_GO(sub_idt_sig$gene)


plt.list = list()
for (i in 1:nrow(sub_idt_sig)){
  plt_df = cbind(t(assay(dt[[1]][sub_idt_sig$gene[i],idx_ftd,'RNA'])), 
                 b_m[idx_ftd,])
  plt_df = as.data.frame(plt_df)
  if (ncol(plt_df)==1){
    print(paste(sub_idt_sig$gene[i], 'missing from the expression set.'))
    next
  }
  colnames(plt_df) = c('Gene', 'Mut')
  plt_df$Gene = as.numeric(plt_df$Gene)
  plt_df_gn = cbind(dt[[1]]@colData[idx_ftd,], plt_df)
  plt_df_gn = as.data.frame(plt_df_gn)
  plt_df_gn$cell_nm = sapply(rownames(plt_df_gn), function(x)return(strsplit(x, split='_')[[1]][1]))
  g = ggplot(plt_df_gn, aes(x=Mut, y=Gene)) + theme_classic() +
    geom_violin() +
    stat_compare_means(comparisons = list(c('mut-','mut+')), method = 't.test') +
    geom_jitter(size=0.5, width=0.2) +
    geom_boxplot(width=0.2, outlier.shape = NA) +
    labs(x=check_mut, y=sub_idt_sig$gene[i]) +
    geom_text(aes(label=cell_nm, alpha=plt_df_gn$Mut=='mut+'), size=2, nudge_x=0, nudge_y=0.02, angle=30) +
    scale_alpha_manual(values = c(0,1)) + 
    theme(legend.position = 'none')
  plt.list[[sub_idt_sig$gene[i]]] = g
}
marrangeGrob(grobs=plt.list,ncol=3,nrow=3) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/eQTL', 'CCLE_hot_spot_from_TCGA_BRCA_dual_sig.pdf'),
         plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)

