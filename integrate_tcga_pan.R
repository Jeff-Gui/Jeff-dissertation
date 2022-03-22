library(tidyverse)
library(maftools)
library(MultiAssayExperiment)
setwd('/Users/jefft/Desktop/p53_project')
source('scripts/eQTL/utils.R')

use_cache = TRUE
if (!use_cache){
  merge_mae = function(mae_a, mae_b){
    rna_a = assay(mae_a[,,'RNA'])
    rna_b = assay(mae_b[,,'RNA'])
    co_gene = intersect(rownames(rna_a), rownames(rna_b))
    rna_b = cbind(rna_a[co_gene,], rna_b[co_gene,])
    meta_merged = rbind(mae_a@colData, mae_b@colData)
    complete_cases = colnames(rna_b)
    exprdat = SummarizedExperiment(rna_b)
    colnames(exprdat) = complete_cases
    multiAssay = MultiAssayExperiment(
      list('RNA'=exprdat),
      meta_merged)
    return(multiAssay)
  }
  
  fle = list.files('datasets')
  coll_rna = list()
  coll_maf = list()
  for (i in fle[grep('[0-9]+-[A-Za-z]+-TCGA', fle)]){
    load(file.path('datasets', i, 'clean_data.RData'))
    ctype = strsplit(i, split='-')[[1]][2]
    sub_maf = dt[[2]]@data
    sub_maf[['cancer_type']] = ctype
    # sub_maf = subset(sub_maf, sub_maf$Hugo_Symbol == 'TP53')
    coll_maf[[ctype]] = sub_maf
    dt[[1]]@colData[['ctype']] = ctype
    coll_rna[[ctype]] = dt[[1]]
  }
  gc()
  coll_maf = merge_mafs(coll_maf)
  merged_mae = Reduce(merge_mae, coll_rna)
  dt = list(merged_mae, coll_maf)
  a = ls()
  rm(list=a[which(a != 'dt')])
  gc()
  
  # leave high tumor grade only
  idx = which(dt[[1]]@colData$AJCC_PATHOLOGIC_TUMOR_STAGE %in% 
    c('STAGE IV','STAGE IVA', 'STAGE IVB', 'STAGE IVC', 'STAGE X'))
  dt[[1]] = dt[[1]][,idx,]
  dt[[2]] = subsetMaf(dt[[2]], tsb = rownames(dt[[1]]@colData))
  
  # save(dt, file = 'datasets/TCGA-Pan-Nine/pan/clean_data.RData')
  save(dt, file = 'datasets/TCGA-Pan-Nine/clean_data.RData')
} else {
  load('datasets/TCGA-Pan-Nine/clean_data.RData')
  # load('datasets/TCGA-Pan-Nine/pan/clean_data.RData')
}

use_cache_PCA = TRUE
library(FactoMineR)
library(factoextra)

if (use_cache_PCA){
  load('datasets/TCGA-Pan-Nine/TCGA_pan_nine_clean_data_pca50.RData')
  ind = get_pca_ind(pca.res)
} else {
  mtx = t(assay(dt[[1]][,,'RNA']))
  pca.res = PCA(mtx, ncp=50, graph = F)
  ind = get_pca_ind(pca.res)
  title = strsplit(config_name, split='\\.')[[1]][1]
  save(pca.res, file = 'datasets/TCGA-Pan-Nine/TCGA_pan_nine_clean_data_pca50.RData')
}
fviz_eig(pca.res, addlabels = T)

library(umap)
library(Rtsne)

# 1. TSNE
pca.tsne = Rtsne(ind$coord, pca = FALSE)
df = pca.tsne$Y
rownames(df) = rownames(ind$coord)

# 2. UMAP
pca.umap = umap(ind$coord)
df = pca.umap
colnames(df) = c('PC1', 'PC2')

# 3. Raw PCA
df = data.frame('PC1'=ind$coord[,1], 'PC2'=ind$coord[,2])

meta = data.frame(dt[[1]]@colData)
df = cbind(df, meta)

# annotate p53
ann = annotate_sample_mut(dt[[2]]@data)
df['p53_state'] = NA
for (i in names(ann)){
  df[ann[[i]], 'p53_state'] = i
}

ggplot(df, aes(x=PC1, y=PC2)) + theme_classic() +
  geom_point(size=1, alpha=0.7, aes(color=as.factor(p53_state))) +
  labs(x='PC1', y='PC2') +
  facet_wrap(~ctype)

# Does high grade tumor has higher p53 mutation rate?
## annotate p53
ann = annotate_sample_mut(dt[[2]]@data)
dt[[1]]@colData['p53_state'] = 'wildtype'
for (i in names(ann)){
  dt[[1]]@colData[ann[[i]], 'p53_state'] = i
}
table(dt[[1]]@colData$p53_state)

meta = as.data.frame(dt[[1]]@colData)
meta = subset(meta, meta$AJCC_PATHOLOGIC_TUMOR_STAGE!='')
map = table(meta$AJCC_PATHOLOGIC_TUMOR_STAGE)
gp = c(1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4)
names(gp) = names(map)
meta$ajcc_stage_numeric = gp[meta$AJCC_PATHOLOGIC_TUMOR_STAGE]
coll = list()
for (i in unique(meta$p53_state)){
  coll[[i]] = meta %>% group_by(ajcc_stage_numeric) %>% 
    summarise(freq=sum(p53_state==i)/length(p53_state)) %>%
    mutate('p53_state'=i)
}
coll = Reduce(rbind, coll)
ggplot(coll, aes(x=ajcc_stage_numeric, y=freq, fill=p53_state)) +
  geom_bar(position = 'stack', stat = 'identity') +
  labs(x='AJCC Tumor Stage', y='Frequency') +
  ggpubr::theme_pubr()
table(meta$ajcc_stage_numeric)

# How many sample lose in Toil?
toil_info = read.table('/Users/jefft/Desktop/p53_project/archived/Zena/TcgaTargetGTEX_phenotype.txt', 
                       sep='\t', header = T)
loss = subset(meta, !rownames(meta) %in% toil_info$sample)
table(loss$ctype) # 171 COAD loss, other almost no loss

toil_sample = read.table('/Users/jefft/Desktop/p53_project/archived/Zena/toil_all_samples.txt')
loss = subset(meta, !rownames(meta) %in% toil_sample$V1)
loss = subset(meta, !toil_sample$V1 %in% rownames(meta))
table(loss$ctype)


