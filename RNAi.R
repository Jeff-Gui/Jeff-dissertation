library(tidyverse)
library(data.table)
library(pheatmap)
setwd('/Users/jefft/Desktop/p53_project')
dt = read.csv('CCLE/RawData/RNAi_breast.csv', header = T, row.names = 1)
dt_no_na = na.omit(t(dt))
dim(dt_no_na)
dist = cor(dt_no_na)
hclust(dist)

meta = as.data.frame(fread('CCLE/RawData/sample_info_21Q3.csv', na.strings = ''))
mts = meta$primary_or_metastasis
names(mts) = meta$DepMap_ID
dt_no_na = dt_no_na[,-which(is.na(mts[colnames(dt_no_na)]))]  # Remove no metastasis record cells 
ann = data.frame('mts'=mts[colnames(dt_no_na)])

pheatmap(dist, show_rownames = F, cluster_rows = T, cluster_cols = T, annotation_col = ann)


# See correlation between p53 mutation and PLK3 RNAi dependency
## Clean PLK3 RNAi data
dt = fread('CCLE/RawData/RNAi_D2_combined_gene_dep_scores.csv', header = T)
dt = as.data.frame(dt)
rownames(dt) = gsub(' \\([0-9]+\\)','',dt$V1)
dt = dt[,-1]
colnames(dt) = gsub('_.*', '', colnames(dt))

plk3_rnai = transpose(dt['PLK3',])
rownames(plk3_rnai) = colnames(dt)
## For RNAi cell line, require it to be present in the sample info
valid_cell = rownames(plk3_rnai)[which(rownames(plk3_rnai) %in% meta$stripped_cell_line_name)]
plk3_rnai = data.frame(plk3_rnai[valid_cell, ])
sub_meta = subset(meta, meta$stripped_cell_line_name %in% valid_cell)
id_to_name = sub_meta$DepMap_ID
names(id_to_name) = sub_meta$stripped_cell_line_name
rownames(plk3_rnai) = id_to_name[valid_cell]
colnames(plk3_rnai) = 'plk3_rnai'

## Read p53 mutation
mut_p53 = read.csv('CCLE/ProcessedData/cl_meta_p53mut_merged.csv', header = T)
ann = c()
for (i in rownames(plk3_rnai)){
  sub_mut = subset(mut_p53, mut_p53$depmap_id==i)
  if (nrow(sub_mut)>=2){
    ann = c(ann, NA)
    next
  }
  if (nrow(sub_mut)==0){
    ann = c(ann, 'no variant')
  } else {
    ann = c(ann, sub_mut$Protein.Change)
  }
}
plk3_rnai['ann'] = ann
plk3_rnai = plk3_rnai[-which(is.na(plk3_rnai['ann'])),]
plk3_rnai = plk3_rnai[-which(is.na(plk3_rnai['plk3_rnai'])),]
plk3_rnai = plk3_rnai[-grep('fs', plk3_rnai$ann),]  # Remove frameshift
plk3_rnai = plk3_rnai[-grep('del', plk3_rnai$ann),]
plk3_rnai = plk3_rnai[-grep('\\*', plk3_rnai$ann),]  # Remove deletion
plk3_rnai = plk3_rnai[-grep('_', plk3_rnai$ann),] # Keep missense only
ann_pos = gsub('p\\.', '', plk3_rnai$ann)
ann_pos = gsub('[A-Z]$', '', ann_pos)
plk3_rnai['pos'] = ann_pos
plk3_rnai = plk3_rnai[order(plk3_rnai$pos),]
table(plk3_rnai$pos)
colnames(plk3_rnai)[1] = 'rnai'
smm = plk3_rnai %>% group_by(pos) %>% summarise(mn=mean(rnai), sd=sd(rnai))
ggplot(smm) + 
  geom_point(aes(x=pos, y=mn)) +
  theme(axis.text.x = element_text(angle=90))


### Specific cancer
meta_sub = subset(meta, meta$primary_disease=='Breast Cancer')
plk3_rnai_sub = plk3_rnai[which(rownames(plk3_rnai_sub) %in% meta_sub$DepMap_ID),]
smm = plk3_rnai_sub %>% group_by(pos) %>% summarise(mn=mean(rnai), sd=sd(rnai))
ggplot(smm) + 
  geom_point(aes(x=pos, y=mn)) +
  theme(axis.text.x = element_text(angle=90))
