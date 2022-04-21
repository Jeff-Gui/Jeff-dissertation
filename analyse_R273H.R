library(tidyverse)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'R273H')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../set_theme.R')
source('../overlap_utils.R')
source('../revigo_utils.R')
library(tidyverse)
library(stringr)
library(ggvenn)


# load gene
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
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
df_coll = subset(df_coll, df_coll$protein_change %in% c('p.R273H'))
pos_gene = subset(df_coll, df_coll$beta > 0)

# Gene overlap: LGG X COAD ====
df_coll$cancer = toupper(sapply(df_coll$experiment, function(x){return(strsplit(x, split='_')[[1]][2])}))
## !!! change sign !!!
df_coll_pos = subset(df_coll, df_coll$beta > 0 & df_coll$cancer != 'STAD') # too few in STAD
gene_coll = list()
for (i in unique(df_coll_pos$cancer)){
  gene_coll[[i]] = df_coll_pos$gene[df_coll_pos$cancer==i]
}
ggvenn(gene_coll)
um = gen_upSet_mtx(gene_coll)
ovl_p = run_overlap_test(um, c(T,T), n_loop = 1000)

pool = rownames(mtx)[mtx$COAD==1 | mtx$LGG==1]
count_exp = c(sapply(list(c(T,T),c(T,F),c(F,T)),function(x){return(length(get_idx_condt(um, x)))}),
  length(!pool %in% rownames(um)))
count_exp = matrix(count_exp,nrow=2)
print(count_exp)
chisq.test(count_exp)

# intersection
check = rownames(um)[get_idx_condt(um, c(T,T))]
test = do_GO(check, background = rownames(um), ont='BP')
dotplot(test)
pcs_GO_out(test, dir = plot_out, filename = 'R273H_BRCA-COAD_neg_co_gene_GO_BP.pdf')

# Cancer specific signature
ont = 'BP'
coll = list()
go_list = list()
for (i in unique(df_coll_pos$cancer)){
  gs = subset(df_coll_pos, df_coll_pos$cancer==i)$gene
  coll[[i]] = gs
  test = do_GO(gs[which(!gs %in% check)], background = unique(df_coll_pos$gene), ont = ont)
  rvg = run_revigo(test@result[,c('ID', 'pvalue')])
  test@result = subset(test@result, ID %in% rvg$`Term ID`)
  go_list[[i]] = pcs_GO_out(test, dir = plot_out, filename = paste('R273H_', i, '_', ont, '_specific.pdf', sep=''))
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, paste('R273H_BRCA-COAD_pos_co_gene.pdf', sep='')),
             plot=., width=4,height=4,units='in',device='pdf',dpi=300)

