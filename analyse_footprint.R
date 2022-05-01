library(tidyverse)
library(readxl)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
output_dir = '/Users/jefft/Desktop/p53_project/datasets/TCGA-ATAC/plots'

## Footprint reported by TCGA-ATAC ====
depth = read_xlsx('/Users/jefft/Desktop/p53_project/datasets/TCGA-ATAC/aav1898_data_s6.xlsx',
                  sheet = 1, skip = 17, col_names = T)
flank = read_xlsx('/Users/jefft/Desktop/p53_project/datasets/TCGA-ATAC/aav1898_data_s6.xlsx',
              sheet = 2, skip = 17, col_names = T)
motif_meta = as.data.frame(depth[,1:4])
rownames(motif_meta) = motif_meta$TF_Motif_Name
length(intersect(motif_meta$TF_Motif_Name, flank$TF_Motif_Name)) # should be same as nrow
depth = as.data.frame(depth)
rownames(depth) = depth$TF_Motif_Name
depth = depth[,-c(1:4)]
flank = as.data.frame(flank)
rownames(flank) = flank$TF_Motif_Name
flank = flank[,-c(1:4)]
print(dim(flank))
save(depth, flank, motif_meta, file = file.path(output_dir, 'clean.RData'))

## Check distribution ====
# values have been log transformed, minus more, deeper; positive more, heigher accessibility.
# see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6408149/
load(file.path(output_dir, 'clean.RData'))
hist(as.numeric(depth[1,]))
hist(as.numeric(flank[1,]))
summary(rowMeans(depth))
# deeper footprint -> lower flanking accessibility
plot(as.numeric(depth[8,]), as.numeric(flank[8,]))

## Correlate to p53 mutation ====
fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData'
dt = load_clean_data(fp, ann_bin_mut_list = c('sandwich', 'conformation', 'contact'))
all_spl = rownames(dt[[1]]@colData)[dt[[1]]@colData$p53_state %in% c('missense', 'Wildtype')]
all_spl = sapply(all_spl, function(x){
  x = strsplit(x, split='-')[[1]]
  x = paste(x[1:3], collapse = '-')
  return(x)
})
all_atac = colnames(depth)
ggvenn::ggvenn(list('ATAC: BRCA + COAD'=all_atac, 'RNA: BRCA'=all_spl))
co_spl = intersect(all_atac, all_spl)
depth_sub = depth[,co_spl]
flank_sub = flank[,co_spl]

dt[[1]] = subsetByColumn(dt[[1]], names(all_spl)[which(all_spl %in% co_spl)])
muts = apply(as.matrix(dt[[1]]@colData[,grep('has_', colnames(dt[[1]]@colData))]), MARGIN = 2, 
      function(x){return(which(x==1))})
for (i in 1:length(muts)){
  assign(strsplit(names(muts)[i], split='_')[[1]][2], names(muts[[i]]))
}
wildtype = rownames(dt[[1]]@colData)[dt[[1]]@colData$p53_state=='Wildtype']
# 42 wt, 7 sandwich, 4 contact, 4 conformation

library(pheatmap)
ann = data.frame('p53_state'=dt[[1]]@colData[,c('p53_state')], 
                row.names = rownames(dt[[1]]@colData))
ann$itg_state = 'other-missense'
ann[contact, 'itg_state'] = 'contact'
ann[conformation, 'itg_state'] = 'conformation'
ann[sandwich, 'itg_state'] = 'sandwich'
ann[wildtype, 'itg_state'] = 'wildtype'
table(ann$itg_state)
rownames(ann) = sapply(rownames(ann), function(x){
  x = strsplit(x, split='-')[[1]]
  x = paste(x[1:3], collapse = '-')
  return(x)
})
ann = ann[order(ann$itg_state),]

UQ = function(m){
  for (i in 1:ncol(m)){
    m[,i] = m[,i] / quantile(m[,i], 0.75, na.rm=T) * 100
  }
  return(m)
}

motif_OI = motif_meta$TF_Motif_Name[grep('Ets', motif_meta$Family_Name)]
motif_OI = motif_meta$TF_Motif_Name[grep('ETS1', motif_meta$TF_Name)]
pheatmap(UQ(flank_sub)[motif_OI,rownames(ann)], show_rownames = F, show_colnames = F, 
         cluster_rows = T, cluster_cols = F, annotation_col = ann, border_color = F,
         scale = 'none')



