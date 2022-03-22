setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
exp_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg'
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)

load(file.path(exp_home, 'GO_result.RData'))
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

# QC gene ratio > 1%
go_df['gene_ratio_val'] = sapply(go_df$GeneRatio, function(x){eval(parse(text = x))})
go_df = go_df[which(go_df$gene_ratio_val > 0.01),]

### Intersect among cancers
lst = list(ids$`tcga_brca_raw_seq-hot_spot-pos`, 
           ids$`tcga_coad_raw_seq-hot_spot-pos`,
           ids$`tcga_lgg_raw_seq-hot_spot-pos`)
names(lst) = c('BRCA', 'COAD', 'LGG')
g = ggvenn(lst)
g %>% ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/', 'LGG_BRCA_COAD_pos_GO.pdf'),
             plot=., width=5,height=5,units='in',device='pdf',dpi=300)

co_term = Reduce(intersect, list(ids$`tcga_brca_raw_seq-hot_spot-pos`, 
                                 ids$`tcga_coad_raw_seq-hot_spot-pos`,
                                 ids$`tcga_lgg_raw_seq-hot_spot-pos`))

co_term_go = subset(hot_spot_go, hot_spot_go$ID %in% co_term & hot_spot_go$mutation == 'hot_spot')

g = ggplot(co_term_go) + theme_classic() +
  geom_point(aes(x=cancer, y=Description, color=-log10(p.adjust), size=Count)) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  mytme +
  labs(title = str_wrap('GO terms of genes positively correlated to p53 hot-spot mutants co-enriched in BRCA, COAD and LGG',
                        width=50))
g %>% ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/', 'GO_co_enrich_hot_spot_BRCA_COAD_LGG.pdf'),
             plot=., width=15,height=8,units='in',device='pdf',dpi=300)

gene_pool = unique(Reduce(c, sapply(co_term_go$geneID, function(x){
  return(strsplit(x, split='/')[[1]])
}))) # a pool of genes belonging to co-enriched terms
output_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg/outputs'
coll = list()
for (i in c('tcga_brca_raw_seq', 'tcga_coad_raw_seq', 'tcga_lgg_raw_seq')){
  fp = file.path(output_home, i, 'trans_eqtl_fdr005.txt')
  dtb = read.table(fp, sep='\t', header = T)
  # dtb['cancer'] = toupper(strsplit(i, split='_')[[1]][2])
  coll[[i]] = dtb$gene
}
coll = Reduce(c, coll)
gene_pool = intersect(coll, gene_pool) # should almost be the same

### Test if they are operational in CCLE
load(file = file.path('/Users/jefft/Desktop/p53_project/datasets/CCLE', 
                      'clean_data_inspect.RData'))
p53_ann = annotate_sample_mut(dt[[2]]@data)
dt[[1]]@colData[['p53_state']] = 'wildtype'
for (i in names(p53_ann)){
  dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
}

out_home = '/Users/jefft/Desktop/p53_project/Plots/Co_enrich_GO'
meta_mut_home = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_mut_fp = list.files(meta_mut_home)
meta_mut_fp = sapply(meta_mut_fp, function(x){return(file.path(meta_mut_home, x))})
meta_mut = load_meta_mut(meta_mut_fp)

idx = rownames(dt[[1]]@colData)
site = sapply(idx, function(x)strsplit(x, split='_')[[1]][2])

##### Set which p53 mutant and which site
check_mut = 'hot_spot'
check_site = c('CENTRAL') # LUNG, BREAST, LARGE, CENTRAL

if (length(grep('p\\.', check_mut))==0){
  mut_aa = unique(meta_mut$aa_pos[meta_mut$meta_mut_id == check_mut])
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
#### Append p53-null information
p53_null_cell = rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state=='nonsense')]
b_m[intersect(rownames(b_m), p53_null_cell),] = 'null'

#### get gene expression
idx_ftd = idx[site %in% check_site] 
gene_pool_expr = t(assay(dt[[1]][gene_pool,idx_ftd,'RNA']))
b_m_ftd = b_m[idx_ftd,]

#### check which specific mutation mut+ cells have
sapply(names(b_m_ftd)[which(b_m_ftd=='mut+')], function(x){
  sub_tb = subset(dt[[2]]@data, dt[[2]]@data$Tumor_Sample_Barcode==x &
                    dt[[2]]@data$Hugo_Symbol=='TP53')
  return(sub_tb$Protein_Change)
})

ps = c()
fcs = c()
nn_count = 0
for (i in 1:length(gene_pool)){
  if (!gene_pool[i] %in% colnames(gene_pool_expr)){
    ps = c(ps, NA)
    fcs = c(fcs, NA)
    print(paste(gene_pool[i], 'not in expression matrix.'))
    nn_count = nn_count + 1
  } else {
    a = gene_pool_expr[b_m_ftd=='mut-',gene_pool[i]]
    b = gene_pool_expr[b_m_ftd=='mut+',gene_pool[i]]
    p = t.test(a,b)$p.value
    ps = c(ps, p)
    fc = mean(b) - mean(a)
    fcs = c(fcs, fc)
  }
}
print(paste(nn_count,'/', length(gene_pool), ' genes not in the expression matrix.', sep=''))
sub_idt = data.frame('gene'=gene_pool, 'ccle.p'=ps, 'ccle.dif'=fcs, 'ccle.FDR' = p.adjust(ps, method = 'fdr')) 
# dual-sig: FDR all <0.05; differene same sign
sub_idt_sig = subset(sub_idt, (ccle.FDR < 0.05) & sub_idt$ccle.dif > 0)

#### GO annotation
# sub_idt_sig[['GO_ann']] = ann_GO(sub_idt_sig$gene)
# sub_idt_sig$GO_ann = as.character(sub_idt_sig$GO_ann)
# write.table(sub_idt_sig, 
#             file = file.path('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg/data_out',
#                              'CCLE_hot_spot_from_TCGA_BRCA_dual_sig.txt'), sep=',',
#             quote = F, row.names = F)


