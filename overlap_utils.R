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
  names(test_condition) = colnames(um)
  test_condition = test_condition[!is.na(test_condition)]
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


load_controls = function(){
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
}




