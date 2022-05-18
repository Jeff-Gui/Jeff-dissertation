library(tidyverse)
setwd('/Users/jefft/Desktop/p53_project')
p53_mut = read.csv('CCLE/cl_meta_p53mut_merged.csv')
# Remove cells with multiple p53 mutations
mpl = p53_mut$depmap_id[duplicated(p53_mut$depmap_id)]
p53_mut = subset(p53_mut, !p53_mut$depmap_id %in% mpl)

wt = subset(p53_mut, is.na(p53_mut$Variant.Type))
mis = subset(p53_mut, p53_mut$Variant.Type == 'SNP' & p53_mut$Variant.Classification == 'Missense_Mutation')

# Analyze by nucleotide change
wt['uid'] = NA
mis['uid'] = paste(mis$Reference.Allele,mis$Start.Position,mis$Tumor.Seq.Allele, sep='')
all_snp = unique(mis$uid)
all_snp = na.omit(all_snp)
pt_changes = c()
for(i in all_snp){
  if (!is.na(i)){
    mis[i] = 0
    wt[i] = 0
    mis[which(mis$uid==i), i] = 1
    pt_change = strsplit(mis[which(mis$uid==i), 'Protein.Change'], split='\\.')[[1]][2]
    pt_changes = c(pt_changes, pt_change)
  }
}
names(pt_changes) = all_snp
coll = rbind(mis, wt)

library(pheatmap)
mtx = as.matrix(coll[,34:ncol(coll)])
colnames(mtx) = pt_changes[colnames(mtx)]
pheatmap(mtx, labels_row = '', cluster_rows = T, cluster_cols = F)

# Analyze by protein change
## We remove non-coding SNP
mis = subset(mis, !is.na(mis$Protein.Change))
wt['pos'] = NA
wt['aa_ori'] = NA
wt['aa_mut'] = NA
mis['pos'] = NA
mis['aa_ori'] = NA
mis['aa_mut'] = NA
pt_change = sapply(mis['Protein.Change'], function(x)return(gsub('p\\.', '', x)))
pt_change = as.vector(na.omit(pt_change))
for(i in pt_change){
  aa_ori = substr(i,1,1)
  if (nchar(i)==5){
    pos = substr(i,2,4)
    aa_mut = substr(i,5,5)
  } else {
    pos = substr(i,2,3)
    aa_mut = substr(i,4,4)
  }
  mis[i] = 0
  wt[i] = 0
  mis[grep(i, mis$Protein.Change), i] = 1
  mis[grep(i, mis$Protein.Change), 'aa_ori'] = aa_ori
  mis[grep(i, mis$Protein.Change), 'aa_mut'] = aa_mut
  mis[grep(i, mis$Protein.Change), 'pos'] = pos
}
mis = mis[order(as.numeric(mis$pos)),]
coll = rbind(mis, wt)

library(pheatmap)
mtx = as.matrix(coll[,36:ncol(coll)])
poss = data.frame(as.numeric(sapply(colnames(mtx), function(x){
  if (nchar(x) == 4){
    return(substr(x,2,3))
  } else {
    return(substr(x,2,4))
  }
})))
rownames(poss) = colnames(mtx)
mtx = mtx[,order(poss[,1])]
rownames(mtx) = coll$depmap_id

row_ann = data.frame('lin' = coll$lineage_1)
rownames(row_ann) = coll$depmap_id
pheatmap(mtx, labels_row = '', cluster_rows = F, cluster_cols = F, 
         filename = 'mut_dist.png', width = 30, height = 10,
         annotation_row = row_ann)

# Check coefficient of Variance of each gene expression
library(data.table)
library(Matrix)
library(Rfast)
prepro = function(mtx, index){
  # A function to pre-process the object from fast read function
  rownm = mtx[,get(index)]
  mtx = mtx[,-1]
  colnm = trimws(gsub('\\(.*\\)', '', colnames(mtx))) # remove () stuff in col name
  colnames(mtx) = colnm
  mtx = as.matrix(mtx)
  rownames(mtx) = rownm
  return(mtx)
}

drop_meta = c('cell_line_display_name', 'lineage_1', 'lineage_2', 'lineage_3', 'lineage_4')
tpm = fread('CCLE/RawData/CCLE_expression_21Q3.csv')
tpm = prepro(tpm, index=colnames(tpm)[1])

# Normalize gene expression with 1. copy number; 2. GAPDH
cnv = fread('CCLE/RawData/CCLE_gene_cn_21Q3.csv')
cnv = prepro(cnv, index='V1')

cons_cell = intersect(rownames(cnv), rownames(tpm))
cons_gene = intersect(colnames(cnv), colnames(tpm))
cnv = cnv[cons_cell, cons_gene]
cnv = round(cnv, 3)
tpm = tpm[cons_cell, cons_gene]
tpm_norm = tpm / cnv
tpm_norm = tpm_norm / tpm_norm[,'GAPDH']
sm = colSums(tpm_norm)
tpm_norm = tpm_norm[,!(is.na(sm) | sm==Inf)]  # Remove genes with TPM but 0 CN
#tpm_norm = tpm / tpm[,'GAPDH']

mean_tpm = colMeans(tpm_norm, na.rm = T)  # omit NA cell lines
sd_tpm = sqrt(colVars(tpm_norm, na.rm = T))
cv = sd_tpm / mean_tpm * 100
dtb = data.frame('mean_log2tpm'=mean_tpm, 'sd_log2tpm'=sd_tpm, 'cv'=cv)
dtb_qc = subset(dtb, dtb$mean_log2tpm>quantile(dtb$mean_log2tpm, 0.05))
hist(dtb_qc$cv, breaks=100)

#df = data.frame(sort(cv)], decreasing = T))
#write.table(df, 'high_cv_genes.rnk', col.names = F, quote = F, sep = '\t')
load('CCLE/ProcessedData/RNASeq_norm.RData')
save(tpm_norm, file = 'CCLE/ProcessedData/RNASeq_norm.RData')

plot(mean_tpm, cv)
ggplot(x=mean_tpm, y=cv) +
  geom_point()

