library(tidyverse)
library(utils)
## Get statistics of gene/GO overlap

gen_upSet_mtx = function(mylist){
  # Generate upSet matrix, mylist should be named with group ID
  pool = unique(unlist(mylist))
  upset_mtx = matrix(0, nrow=length(pool), ncol=length(mylist))
  rownames(upset_mtx) = pool
  colnames(upset_mtx) = names(mylist)
  upset_mtx = as.data.frame(upset_mtx)
  for (i in names(mylist)){
    hits = mylist[[i]]
    upset_mtx[hits, i] = 1
  }
  to_rm = which(rowSums(upset_mtx)==0)
  if (length(to_rm) > 0){
    upset_mtx = upset_mtx[-to_rm,]
  }
  return(upset_mtx)
}


get_idx_condt = function(um, test_condition){
  if (is.null(names(test_condition))){
    names(test_condition) = colnames(um)
    test_condition = test_condition[!is.na(test_condition)]
  } else {
    tcd = rep(NA, ncol(um))
    names(tcd) = colnames(um)
    for (i in 1:length(test_condition)){
      tcd[names(test_condition)[i]] = test_condition[i]
    }
    test_condition = tcd[!is.na(tcd)]
  }
  idx = 1:nrow(um)
  #### get observed value (overlap size)
  for (j in 1:length(test_condition)){
    tc = test_condition[j]
    tcnm = names(test_condition)[j]
    if (tc){
      idx = intersect(idx, which(um[,tcnm]==1))
    } else {
      idx = intersect(idx, which(um[,tcnm]==0))
    }
  }
  return(idx)
}


run_overlap_test = function(um, test_condition, n_loop=1000){
  # Shuffle upset matrix, how many observed is significant?
  # um: upset binary matrix (feature ~ group)
  # test_condition: same length as the ncol of upset matrix, 
  #   code with the following: T: must be 1, F: must be 0, NA: no matter. 
  pb = txtProgressBar(style=3)
  names(test_condition) = colnames(um)
  idx = get_idx_condt(um, test_condition)
  obs = length(idx)
  test_condition = test_condition[!is.na(test_condition)]
  #### get null distribution: relabel each column without replacement
  null_dist = c()
  for (i in 1:n_loop){
    idx = 1:nrow(um)
    test = apply(um, MARGIN = 2, FUN = function(x){sample(x,replace = F)})
    for (j in 1:length(test_condition)){
      tc = test_condition[j]
      tcnm = names(test_condition)[j]
      if (tc){
        idx = intersect(idx, which(test[,tcnm]==1))
      } else {
        idx = intersect(idx, which(test[,tcnm]==0))
      }
    }
    null_dist = c(null_dist, length(idx))
    setTxtProgressBar(pb, i/n_loop)
  }
  close(pb)
  return(list('observed'=obs, 'null_dst'=null_dist, 'p.left'=sum(obs > null_dist) / length(null_dist)))
}


gen_list_from_um = function(um){
  rt = list()
  for (i in 1:ncol(um)){
    rt[[colnames(um)[i]]] = rownames(um)[um[,i]==1]
  }
  return(rt)
}

append_um_col = function(mtx, gene_list, list_name){
  # append a column to upset matrix
  rt = gen_list_from_um(mtx)
  rt = append(rt, list(gene_list))
  names(rt)[length(names(rt))] = list_name
  return(gen_upSet_mtx(rt))
}

load_controls = function(um_to_merge=NULL, ng_ctrs=FALSE, sel_mut=NULL){
  #### load pre-computed CCLE BRCA test
  ccle.brca = read.table('/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd/BRCA/hotspot_comp.txt', sep='\t', header = T)
  ccle.brca.mask = ccle.brca
  ccle.brca.mask[is.na(ccle.brca.mask)] = 0
  ccle.rna.up = ccle.brca.mask$gene[ccle.brca.mask$rna_hs.wt_dif > 0 & ccle.brca.mask$rna_hs.wt_p < 0.05]
  ccle.rnai.down = ccle.brca.mask$gene[ccle.brca.mask$rnai_hs.wt_dif < 0 & ccle.brca.mask$rnai_hs.wt_p < 0.05]
  ccle.crispr.down = ccle.brca.mask$gene[ccle.brca.mask$crispr_hs.wt_dif < 0 & ccle.brca.mask$crispr_hs.wt_p < 0.05] # nothing if use padj
  ccle.rna.up.null = ccle.brca.mask$gene[ccle.brca.mask$rna_hs.ns_dif > 0 & ccle.brca.mask$rna_hs.ns_p < 0.05]
  #### load literature control
  known_target = read.table('/Users/jefft/Desktop/p53_project/datasets/Fischer2017/TableS2.txt', header = T)
  known_wt_up = known_target$Gene.Symbol[known_target$sum.direct.regulation.score...16..0..16. >= 2]
  known_wt_down = known_target$Gene.Symbol[known_target$sum.direct.regulation.score...16..0..16. <= -2]
  #### load H1299 R273H ChIP-seq
  peaks = read.table('/Users/jefft/Desktop/p53_project/datasets/PRJEB20314/p53-1-2-merged_Peaks.txt', header = T)
  mean_chek1 = mean(peaks$mean_ern.rep1[peaks$gene=='CHEK1'], peaks$mean_ern.rep2[peaks$gene=='CHEK1'])
  peaks = na.omit(peaks)
  peak_over_chek1 = peaks$gene[rowMeans(as.matrix(peaks[,c('mean_ern.rep1', 'mean_ern.rep2')])) >= mean_chek1]
  
  rt_list = list('ccle.rna.up'=ccle.rna.up, 'ccle.rnai.down'=ccle.rnai.down, 'ccle.rnai.up.null'=ccle.rna.up.null,
                 'ccle.crispr.down' = ccle.crispr.down,
                 'known.wt.up'=known_wt_up, 'known.wt.down'=known_wt_down, 'peak.over.chek1'=peak_over_chek1)
  
  #### load negative controls
  if (ng_ctrs){
    library(readxl)
    ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
    ctrs = na.omit(ctrs_raw)
    if (!is.null(sel_mut)){
      ctrs = ctrs[ctrs$Gene %in% sel_mut]
    }
    ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
    ctrs = ctrs[-which(duplicated(ctrs)),]
    ng1 = ctrs$Gene[ctrs$Gene_annotation=='wt_control']
    ng2 = ctrs$Gene[ctrs$Gene_annotation=='wt_control_2']
    pg = ctrs$Gene[ctrs$Gene_annotation!='wt_control' & ctrs$Gene_annotation!='wt_control_2']
    if (length(ng1)>0){
      rt_list[['ng1']] = ng1
    }
    if (length(ng2)>0){
      rt_list[['ng2']] = ng2
    }
    if (length(pg)>0){
      rt_list[['pg']] = pg
    }
  }
  
  if (!is.null(um_to_merge)){
    rt_list = append(rt_list, gen_list_from_um(um_to_merge))
    return(gen_upSet_mtx(rt_list))
  } else {
    return(rt_list)
  }
}




