## Do hypothesis testing in CCLE data
library(tidyverse)
library(ggpubr)
library(utils)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('enrich_utils.R')
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')


# Load data ====
CCLE_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
load(file.path(CCLE_home, 'BRCA', 'clean_data.RData'))
rnai = load_rnai()
crispr = load_crispr_score()

# Compare hotspot to wildtype, RNA/RNAi/CRISPR ====
mut_ann = annotate_mut_group(dt, check_muts = c('hot_spot', 'conformation', 'contact', 'sandwich'), 
                   as.mtx = TRUE)
hotspot = rownames(mut_ann)[mut_ann$hot_spot==1]
wildtype = rownames(mut_ann)[mut_ann$wildtype==1]
nonsense = rownames(mut_ann)[mut_ann$nonsense==1]

hs_rnai = intersect(colnames(rnai), hotspot)
hs_crispr = intersect(colnames(crispr), hotspot)
wt_rnai = intersect(colnames(rnai), wildtype)
wt_crispr = intersect(colnames(crispr), wildtype)

# Run testing ====
rna = assay(dt[[1]][,,'RNA'])

t.test.wrap = function(a,b){
  a = na.omit(as.numeric(a))
  b = na.omit(as.numeric(b))
  dif = mean(a) - mean(b)
  p = t.test(a,b)$p.value
  return(list(dif, p))
}

run_bulk_t = function(df, group1, group2, p_mtd='fdr'){
  pb = txtProgressBar(style=3)
  v = vector(length=nrow(df))
  rsl = list('gene'=v,'dif'=v, 'p'=v)
  for (i in 1:nrow(df)){
    gene = rownames(df)[i]
    rs = t.test.wrap(df[gene, group1], df[gene, group2])
    rsl[['gene']][i] = gene
    rsl[['dif']][i] = rs[[1]]
    rsl[['p']][i] = rs[[2]]
    setTxtProgressBar(pb, i/nrow(df))
  }
  rsl[['p.adj']] = p.adjust(rsl[['p']], method = p_mtd)
  rsl = as.data.frame(rsl)
  close(pb)
  return(rsl)
}

result = list()
### RNA: Hotspot VS Wildtype ====
result[['rna_hs.wt']] = run_bulk_t(rna, hotspot, wildtype)

### RNA: Hotspot VS Nonsense ====
result[['rna_hs.ns']] = run_bulk_t(rna, hotspot, nonsense)

### RNAi: Hotspot VS Wildtype ====
result[['rnai_hs.wt']] = run_bulk_t(rnai, hs_rnai, wt_rnai)

### CRISPR: Hotspot VS Wildtype ====
result[['crispr_hs.wt']] = run_bulk_t(crispr, hs_crispr, wt_crispr)

## Add dummy NA rows
pool = sort(unique(c(rownames(rnai), rownames(crispr), rownames(rna))))
for (i in names(result)){
  to_add = pool[which(!pool %in% result[[i]]$gene)]
  if (length(to_add) > 0){
    addm = matrix(NA, nrow=length(to_add), ncol=ncol(result[[i]]))
    addm = as.data.frame(addm)
    colnames(addm) = colnames(result[[i]])
    addm$gene = to_add
    result[[i]] = rbind(result[[i]], addm)
  }
  rownames(result[[i]]) = result[[i]]$gene
  result[[i]] = result[[i]][,-which(colnames(result[[i]])=='gene')]
  result[[i]] = result[[i]][pool,]
  colnames(result[[i]]) = paste(i, colnames(result[[i]]), sep='_')
}
result_coll = Reduce(cbind, result)
result_coll$gene = rownames(result_coll)
write.table(result_coll, file.path(CCLE_home, 'BRCA/hotspot_comp.txt'), 
            sep='\t', row.names = F, quote = F, na = '')


# Archived ====
union_gene = unique(c(rownames(crispr), rownames(rnai), rownames(crispr)))
result = matrix(NA, nrow=length(union_gene), ncol=3*3)
colnames(result) = c('rna.dif', 'rna.p', 'rna.adj.p',
                     'rnai.dif', 'rnai.p', 'rnai.adj.p',
                     'crispr.dif', 'crispr.p', 'crispr.adj.p')
rownames(result) = union_gene
result = as.data.frame(result)

for (i in 1:nrow(result)){
  gene = rownames(result)[i]
  if (gene %in% rownames(rna)){
    rs = t.test.wrap(rna[gene, hotspot], rna[gene, wildtype])
    result[i, 'rna.dif'] = rs[[1]]
    result[i, 'rna.p'] = rs[[2]]
  }
  if (gene %in% rownames(rnai)){
    rs = t.test.wrap(rnai[gene, hs_rnai], rnai[gene, wt_rnai])
    result[i, 'rnai.dif'] = rs[[1]]
    result[i, 'rnai.p'] = rs[[2]]
  }
  if (gene %in% rownames(crispr)){
    rs = t.test.wrap(crispr[gene, hs_crispr], crispr[gene, wt_crispr])
    result[i, 'crispr.dif'] = rs[[1]]
    result[i, 'crispr.p'] = rs[[2]]
  }
  
  if (i %% 2000 == 0){
    print(i)
  }
}
gc()
result$gene = rownames(result)
write.table(result, file.path(CCLE_home, 'BRCA/hotspot_VS_wt.txt'), 
            sep='\t', row.names = F, quote = F, na = '')

### Only keep consistent genes
triple_full = intersect(which(!is.na(result$rnai.dif)), which(!is.na(result$crispr.dif)))
triple_full = intersect(triple_full, which(!is.na(result$rna.dif)))
result = result[triple_full,]

result = result[which((result$rnai.dif * result$crispr.dif > 0) & 
                      (result$rnai.dif * result$rna.dif < 0)),]

for (i in c(3,6,9)){
  adjp = p.adjust(na.omit(result[,i-1]), method = 'fdr')
  j = 1
  t = 1
  while (j <= nrow(result)){
    if (!is.na(result[j,i-1])){
      result[j,i] = adjp[t]
      t = t +1
    }
    j = j + 1
  }
}

result = result[order(result$rnai.dif),]
