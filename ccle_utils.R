library(tidyverse)
library(data.table)

load_rnai = function(rnai_fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE/extra_raw/RNAi_D2_combined_gene_dep_scores.csv',
                     hg19_ann = '/Users/jefft/Desktop/p53_project/scripts/eQTL/hg19_gene_table_autosome.tsv',
                     protein_coding_only = TRUE){
  tb = fread(rnai_fp, sep=',', header = T, na.strings = 'NA', quote = '\"')
  tb = as.data.frame(tb)
  rownames(tb) = sapply(tb$V1, function(x){
    return(strsplit(x, split = ' \\(')[[1]][1])
  })
  tb = tb[,-1]
  if (protein_coding_only){
    # CCLE is aligned to hg19
    genepos = read.table(hg19_ann, sep='\t', header = T)
    tb = tb[which(rownames(tb) %in% genepos$geneid),]
  }
  tb = tb[order(rownames(tb)),]
  return(tb)
}


load_protein_z_score = function(protein_fp='/Users/jefft/Desktop/p53_project/datasets/CCLE/ccle_broad_2019/data_protein_quantification_zscores.txt'){
  tb = fread(protein_fp, sep='\t', header = T, na.strings = 'NA', quote = '\"')
  tb = as.data.frame(tb)
  rownames(tb) = tb[,1]
  tb = tb[,-1]
  # collapse duplicated proteins to per-gene by taking mean
  gene_nms = sapply(rownames(tb), function(x){
    return(strsplit(x, split = '\\|')[[1]][1])
  })
  tb['gene_names'] = gene_nms
  dup_idx = which(duplicated(gene_nms))
  dup_nm = gene_nms[dup_idx]
  tb_non_dup = tb[-which(tb$gene_name %in% dup_nm),]
  tb_dup = tb[which(tb$gene_name %in% dup_nm),]
  tb_dedup = gather(tb_dup, key='cl', value='value', 1:(ncol(tb_dup)-1)) %>%
    group_by(gene_names, cl) %>% summarise(mean_value = mean(value, na.rm=T)) %>%
    spread(cl, mean_value)
  tb_dedup = as.data.frame(tb_dedup)
  rownames(tb_dedup) = tb_dedup$gene_names
  tb_dedup = tb_dedup[,-1]
  tb_dedup[is.nan(as.matrix(tb_dedup))] = NA
  tb_dedup['gene_names'] = rownames(tb_dedup)
  tb_dedup = tb_dedup[,colnames(tb)]
  
  # re-calculate the z-score before merging?
  tb_non_dup = rbind(tb_non_dup, tb_dedup)
  tb_non_dup = tb_non_dup[order(rownames(tb_non_dup)),]
  rownames(tb_non_dup) = tb_non_dup$gene_names
  tb_non_dup = tb_non_dup[,-which(colnames(tb_non_dup)=='gene_names')]
  return(tb_non_dup)
}

