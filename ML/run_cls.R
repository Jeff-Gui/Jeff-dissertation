# TCGA-pan_VS-mutneg_ult, TCGA-pan_VS-wt, TCGA-pan_VS-null
data_out = '/Users/jefft/Desktop/p53_project/ML_experiments/outputs'
setwd('/Users/jefft/Desktop/p53_project/scripts/ML')
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
source('../eQTL/utils.R')
source('../set_theme.R')

source('../eQTL/enrich_utils.R')
setwd('/Users/jefft/Desktop/p53_project/scripts/ML')

# Clean data ====
fp = '/Users/jefft/Desktop/p53_project/datasets/METABRIC/clean_data_new.RData'
fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data_noRankNorm.RData'

## For METABRIC data
# dt = load_clean_data(fp = fp, ann_bin_mut_list = c('hot_spot'), mode='tcga')
# library(data.table)
# source('../eQTL/load_data_cbp.R')
# exp_mtx = fread('/Users/jefft/Desktop/p53_project/datasets/METABRIC/brca_metabric/data_mrna_agilent_microarray_zscores_ref_all_samples.txt')
# exp_mtx = as.data.frame(exp_mtx)
# complete_cases = read.table('/Users/jefft/Desktop/p53_project/datasets/METABRIC/brca_metabric/case_lists/cases_complete.txt', sep=':')
# complete_cases = trimws(strsplit(complete_cases[nrow(complete_cases),2], '\t')[[1]])
# exp_mtx = clean_matrix(exp_mtx, complete_cases_dot = complete_cases, gene_col_nm = 'Hugo_Symbol')
# genepos = read.table('/Users/jefft/Desktop/p53_project/scripts/eQTL/hg18_gene_table_autosome.tsv', header = T)
# exp_mtx = exp_mtx[intersect(rownames(exp_mtx), genepos$geneid),]
# exprdat = SummarizedExperiment(exp_mtx)
# meta = dt[[1]]@colData
# dt[[1]] = MultiAssayExperiment(list('RNA'=exprdat),meta)
# save(dt, file='/Users/jefft/Desktop/p53_project/datasets/METABRIC/clean_data_zScoreRefAll.RData')

# metabric
load('/Users/jefft/Desktop/p53_project/datasets/METABRIC/clean_data_zScoreRefAll.RData')
# tcga
dt = load_clean_data(fp = fp, ann_bin_mut_list = c('hot_spot'), mode='tcga')

idx = which(dt[[1]]@colData$has_hot_spot==1)
idx = c(idx, which(dt[[1]]@colData$p53_state %in% c('nonsense', 'Wildtype')))
idx = which(dt[[1]]@colData$p53_state %in% c('missense', 'nonsense'))
meta = as.data.frame(dt[[1]]@colData[idx,])
meta = meta[!meta$p53_state %in% c('others', 'frameshift'),]
X = t(assay(dt[[1]][,rownames(meta),'RNA']))
nmmap = c(1,2,3)
names(nmmap) = unique(meta$p53_state)
# retrieve VAF as sample weights
maf_idx = which(dt[[2]]@data$Tumor_Sample_Barcode %in% rownames(meta))
maf = as.data.frame(dt[[2]]@data[maf_idx,])
maf = maf[maf$Hugo_Symbol=='TP53',]
rownames(maf) = maf$Tumor_Sample_Barcode
sample_w = maf[rownames(meta),'VAF']
meta$HGVSp_Short = maf[rownames(meta), 'HGVSp_Short']

y = nmmap[meta$p53_state]
print(table(y))

## only train with nonsense and hotspot???
nshs = which(names(y) %in% c('missense','nonsense'))
y = y[nshs]
X = X[nshs,]
nmmap = nmmap[c(1,3)]

names(y) = NULL

# Visualise input ====
library(FactoMineR)
library(factoextra)
X = X[,top2000]
pca.res = PCA(scale(X), 
              ncp=20, graph = F)
ind = get_pca_ind(pca.res)
scre_plot = fviz_eig(pca.res, addlabels = T)
scre_plot
#ggsave(file.path(out_dic, 'genotype_PCA.png'), scre_plot, device = 'png')
#df = data.frame('PC1'=ind$coord[,1], 'PC3'=ind$coord[,20])
#ggplot(df) +
#  geom_point(aes(x=PC1, y=PC3))
pca_coord = as.data.frame(ind$coord)
pca_coord = cbind(pca_coord, dt[[1]]@colData[rownames(pca_coord),])
# pca_coord_sub = pca_coord[pca_coord$ER_STATUS=='Positive'&pca_coord$HER2_STATUS=='Negative',]
ggplot(pca_coord, aes(x=Dim.3, y=Dim.4)) +
  geom_point(size=0.5, aes(color=p53_state)) +
  mytme

df_exp = data.frame('state'=meta$p53_state, X[,top2000])
ps = c()
dif = c()
for (i in top2000){
  a = df_exp[df_exp$state=='missense',i]
  b = df_exp[df_exp$state=='nonsense',i]
  ps = c(ps, wilcox.test(a, b)$p.value)
  dif = c(dif, mean(a) - mean(b))
}
deg = data.frame('gene'=top2000, 'p.value'=p.adjust(ps, method = 'fdr'), 'dif'=dif)

ggplot(gather(df_exp[,c('state','TP53')], key='gene',value='vl',2), aes(x=state,y=vl)) +
  geom_violin() +
  stat_compare_means() +
  facet_wrap(~gene)

# Clustering ====
fviz_nbclust(X, kmeans, method='wss', k.max = 20)
hc = hclust(d = dist(X, method = 'euclidean'), method = 'ward.D2')
fviz_dend(hc, cex=0.2, k=3)
cst_id = cutree(hc, k=3)
meta$cst_id = cst_id[rownames(meta)]
table(meta[meta$p53_state=='nonsense','cst_id'])
table(meta[meta$p53_state=='missense','cst_id'])
table(meta[meta$has_hot_spot==1,'cst_id'])
table(meta[meta$cst_id==3,'HGVSp_Short'])

# Prepare Python env ====
library(reticulate)
use_condaenv('/opt/anaconda3/envs/py37')
source_python('ml_model.py')

# Optimize a random forest model ====
test = run_opt(X,y,cv_round=5,sample_weight=sample_w) # 180 fit, 2.9 min on mac (TCGA), 8.5 min on mac (METABRIC)
#saveRDS(X, file='X.rds')
#saveRDS(y, file='y.rds')
#   param_grid={'max_depth': array([ 4,  6,  8, 10, 12]),
# 'min_samples_leaf': [4, 8, 16, 32],
# 'n_estimators': array([100, 400, 700])},
# Best: {'max_depth': 8, 'min_samples_leaf': 4, 'n_estimators': 100}
# use roc to evaluate
# Best: {'max_depth': 10, 'min_samples_leaf': 16, 'n_estimators': 100}
print(test)

model_suffix = 'METABRIC-Dual'
model_suffix = 'TCGA-NSMS-SVM'

# forest
cls_rs = run_fit_forest(X,y,list('n_estimators'=100, 'max_depth'=6, 'min_samples_leaf'=4),
                 save_fp = file.path(data_out, paste('external_',model_suffix,'.pkl', sep='')))
print(cls_rs$oob)
# SVM
cls_rs = run_fit_SVM(X,y,list('C'=0.001),
                 save_fp = file.path(data_out, paste('external_',model_suffix,'.pkl', sep='')))

inverse_map = function(x){
  # inverse the name and value of a vector
  y = names(x)
  names(y) = x
  return(y)
}

get_rf_plt = function(cls_rs){
  pred_gp = apply(cls_rs$predict_mtx, MARGIN = 1, function(x){return(which(x==max(x)))})
  pred_score = apply(cls_rs$predict_mtx, MARGIN = 1, function(x){return(max(x))})
  out = data.frame('gt'=y, 'pred'=pred_gp, 'pred_score' = cls_rs$predict_mtx, 'pred_max_score'=pred_score)
  out$sample = rownames(X)
  out = out[order(out$gt),]
  inv_nmmap = inverse_map(nmmap)
  out$gt = inv_nmmap[out$gt]
  out$pred = inv_nmmap[out$pred]
  return(out)
}

myPalette = colorRampPalette(c("royalblue","purple", "coral"))
sc = scale_colour_gradientn(colours = myPalette(100), name='Prediction probability', 
                            breaks = c(0.4,0.7,1), limits=c(0,1))

plt_df = get_rf_plt(cls_rs)
xlab = unique(plt_df$pred)
g = ggplot(plt_df, aes(x=gt, y=pred)) +
  geom_jitter(size=0.5, alpha=0.7, aes(color=pred_max_score), width = 0.3) +
  labs(y='Predicted',x='Ground truth') +
  geom_vline(xintercept = c(1.5,2.5), color='grey') +
  geom_hline(yintercept = c(1.5,2.5), color='grey') +
  scale_x_discrete(label = xlab) +
  scale_y_discrete(label = xlab) +
  sc + mytme +
  theme(legend.text = element_text(size=10, face='bold'), legend.direction = 'horizontal')
g
ggsave(file.path(data_out, 'external_cls_GT.pdf'),bg='transparent',
       plot=g, height=11.69*0.4,width=8.27*0.6,units='in',device='pdf',dpi=300)

# Generate the secondary model ====
fimp = t(cls_rs$coef)
fimp = cls_rs$feature_importance
hist(fimp,breaks=100)
quantile(fimp, 0.95)
search_range = seq(10,2500,10)
cum_imp = c()
for (i in search_range){
  cum_imp = c(cum_imp, sum(sort(fimp,decreasing = T)[1:i]))
}
top_n = 900 # feature size of the internal classifier
g = ggplot(data.frame('top_n'=search_range,'prob_cml'=cum_imp)) +
  mytme +
  geom_vline(xintercept = top_n, color='red', size=0.5) +
  geom_hline(yintercept = sum(sort(fimp,decreasing = T)[1:top_n]), color='red',size=0.5) +
  geom_point(aes(x=top_n, y=prob_cml), size=0.5) +
  scale_y_continuous(breaks = seq(0.1,1,0.1)) +
  labs(x='Top n genes', y='Cumulative importance')
g
ggsave(file.path(data_out, 'external_feature_importance.pdf'),bg='transparent',
       plot=g, height=11.69*0.4,width=8.27*0.6,units='in',device='pdf',dpi=300)

print(sum(sort(fimp,decreasing = T)[1:900]))

top2000 = order(fimp, decreasing = T)[1:top_n]
top2000 = colnames(X)[top2000]
# enrichment high important genes
test = do_GO(top2000[1:200])

# saveRDS(X[,top2000], file='X_top2k.rds')
test = run_opt(X[,top2000],y,cv_round=5) # 180 fit, 3 min on mac
print(test)
cls_rs_2 = run_fit_forest(X[,top2000],y,list('n_estimators'=100, 'max_depth'=8, 'min_samples_leaf'=4),
                   save_fp = file.path(data_out, paste('internal_',model_suffix,'.pkl', sep='')))
cls_rs_2 = run_fit_forest(X[,top2000],y,list('n_estimators'=500, 'max_depth'=8, 'min_samples_leaf'=4),
                   save_fp = file.path(data_out, paste('internal_',model_suffix,'.pkl', sep='')))
print(cls_rs_2$oob)
plt_df = get_rf_plt(cls_rs_2)
g = ggplot(plt_df, aes(x=gt, y=pred)) +
  geom_jitter(size=0.5, alpha=0.7, aes(color=pred_max_score), width = 0.3) +
  labs(y='Predicted',x='Ground truth') +
  geom_vline(xintercept = c(1.5,2.5), color='grey') +
  geom_hline(yintercept = c(1.5,2.5), color='grey') +
  scale_x_discrete(label = xlab) +
  scale_y_discrete(label = xlab) +
  sc + mytme +
  theme(legend.text = element_text(size=10, face='bold'), legend.direction = 'horizontal')
g
ggsave(file.path(data_out, 'internal_cls_GT.pdf'),bg='transparent',
       plot=g, height=11.69*0.4,width=8.27*0.6,units='in',device='pdf',dpi=300)

# Run in explorative set ====
X_sample = setdiff(rownames(dt[[1]]@colData)[dt[[1]]@colData$p53_state=='missense'],
                   rownames(X))
X_test = assay(dt[[1]][top2000,X_sample,'RNA'])
# X_test = assay(dt[[1]][top2000,,'RNA'])  # spike in all training
# cls_rs_test = run_cls(t(assay(dt[[1]][,X_sample,'RNA'])),
#                      model_fp = file.path(data_out, paste('internal_',model_suffix,'.pkl', sep='')))
cls_rs_test = run_cls(t(X_test), 
                      model_fp = file.path(data_out, paste('internal_',model_suffix,'.pkl', sep='')))

pred_gp = apply(cls_rs_test, MARGIN = 1, function(x){return(which(x==max(x)))})
pred_score = apply(cls_rs_test, MARGIN = 1, function(x){return(max(x))})
inv_nmmap = inverse_map(nmmap)
cls_test_df = data.frame(cls_rs_test, 'group'=inv_nmmap[pred_gp], 'top_score' = pred_score)
# cls_test_df$sample = rownames(dt[[1]]@colData)
cls_test_df$sample = X_sample
cls_test_df$group = gsub('Wildtype', 'WT-like', cls_test_df$group)
cls_test_df$group = gsub('missense', 'HS-like', cls_test_df$group)
cls_test_df$group = gsub('nonsense', 'NS-like', cls_test_df$group)
table(cls_test_df$group)

ggplot(cls_test_df %>% group_by(group) %>% summarise(count=length(top_score))) +
  geom_histogram(aes(x=1,fill=group,y=count), stat = 'identity') +
  coord_polar(theta = 'y') +
  mytme

save(cls_test_df, cls_rs, cls_rs_2, meta, file = file.path(data_out, 'out.RData'))

# Analyse correlation with MAF SIFT prediction ====
load('/Users/jefft/Desktop/p53_project/datasets/UMD_TP53/SIFT_UMD.RData')
subMaf = subsetMaf(dt[[2]], tsb = cls_test_df$sample, genes = 'TP53')
subMaf = subMaf@data
subMaf$SIFT_UMD = sift_pred[subMaf$HGVSp_Short]
table(subMaf$SIFT_UMD)

# Analyse correlation with clinical features ====
library(ggpubr)
meta = as.data.frame(dt[[1]]@colData)
meta_X = meta
meta_X$p53_state[meta_X$has_hot_spot==1 & meta_X$p53_state=='missense'] = 'hotspot'
## p53 mutation correlated to aneuploidy score, known, but no difference among hotspot, missense and nonsense.
ggplot(meta_X, aes(x=as.factor(p53_state), y=ANEUPLOIDY_SCORE)) +
  geom_violin() + geom_boxplot(width=0.2) +
  stat_compare_means()

cls_test_df_sub = cls_test_df[cls_test_df$top_score > 0.5,]
meta = meta[cls_test_df_sub$sample,]
meta[cls_test_df_sub$sample,'predict_group'] = cls_test_df_sub$group
num_var_to_test = c()
ggplot(meta, aes(x=predict_group, y=ANEUPLOIDY_SCORE)) +
  geom_violin() + geom_boxplot(width=0.2) + stat_compare_means()


# Analyse the maf / mutation status of the explorative set ====
maf_expl = subsetMaf(dt[[2]], genes = 'TP53', tsb = cls_test_df$sample)
maf_expl = maf_expl@data
maf_wt_like = maf_expl[maf_expl$Tumor_Sample_Barcode %in% cls_test_df$sample[cls_test_df$group=='WT-like']]
maf_ns_like = maf_expl[maf_expl$Tumor_Sample_Barcode %in% cls_test_df$sample[cls_test_df$group=='NS-like']]
maf_hs_like = maf_expl[maf_expl$Tumor_Sample_Barcode %in% cls_test_df$sample[cls_test_df$group=='HS-like']]
freq = list()

for (i in unique(cls_test_df$group)){
  msub = maf_expl[maf_expl$Tumor_Sample_Barcode %in% cls_test_df$sample[cls_test_df$group==i]] 
  msub = as.data.frame(table(msub$Protein_position))
  colnames(msub) = c('pos', paste('count', strsplit(i, split='-')[[1]][1], sep='_'))
  freq[[i]] = msub
}
freq_joined = Reduce(full_join, freq)
freq_joined[is.na(freq_joined)] = 0
for (i in 2:ncol(freq_joined)){
  freq_joined[,i] = as.numeric(freq_joined[,i])
}
#freq_joined = freq_joined[which(rowSums(freq_joined[,2:4])>5),]
library(ggrepel)
g = ggplot(freq_joined, aes(x=freq_joined$count_HS, y=freq_joined$count_WT)) +
  geom_jitter(width = 0.01) +
  geom_abline(intercept = 0, slope=1) + mytme +
  geom_text_repel(aes(label=pos))
g
ggsave(file.path(data_out, 'single_mut_prediction.pdf'),bg='transparent',
       plot=g, height=11.69*0.4,width=8.27*0.6,units='in',device='pdf',dpi=300)

meta = as.data.frame(dt[[1]]@colData)
cls_test_df_sorted = cls_test_df[order(cls_test_df$group),]
vis_exp = cbind(X_test[1:100,cls_test_df_sorted$sample], 
                (t(as.data.frame(X))[rownames(X_test)[1:100],]))
ann_col = data.frame('group' = c(cls_test_df_sorted$group, inv_nmmap[y]), 
                     row.names = colnames(vis_exp),
                     'subtype'=meta[colnames(vis_exp),'CLAUDIN_SUBTYPE'])
ann_col = ann_col[order(ann_col$group, ann_col$subtype),]
ann_col$group = factor(ann_col$group, levels = c('Wildtype', 'WT-like', 'missense', 'HS-like', 'nonsense', 'NS-like'))
mn = floor(min(vis_exp))
mx = ceiling(max(vis_exp))
bk = c(seq(mn,-0.1,by=0.1),seq(0,mx,by=0.1))
pheatmap(vis_exp[rownames(ann_col)], cluster_rows = T, cluster_cols = F, 
         show_rownames = F, show_colnames = F, scale='none',
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks = bk, legend_breaks = seq(mn,mx,2),
         annotation_col = ann_col, filename = file.path(data_out, 'heatmap.pdf'))



