library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(gridExtra)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')

ann_GO = function(gene_list, category='BP'){
  out_list = list()
  gene_ann = AnnotationDbi::select(org.Hs.eg.db,
                    columns = c('SYMBOL', 'GO'),
                    keys = keys(org.Hs.eg.db, keytype = 'ENTREZID'))
  gene_ann_ftd = subset(gene_ann, gene_ann$ONTOLOGY == category &
                        gene_ann$SYMBOL %in% gene_list)
  # filter experimental evidence
  gene_ann_ftd = subset(gene_ann_ftd, gene_ann_ftd$EVIDENCE %in% c(
    'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'
  ))
  for (i in gene_list){
    terms = unlist(subset(gene_ann_ftd, gene_ann_ftd$SYMBOL==i)$GO)
    out_list[[i]] = paste(sapply(terms, function(x)return(GOTERM[[x]]@Term)), collapse = '; ')
  }
  return(out_list)
}

do_GO = function(df){
  avg_score = aggregate(beta~gene, df, mean)
  avg_score = avg_score[order(avg_score$beta, decreasing = T),]
  overlapping_gene = avg_score$gene
  gene_ETR = bitr(overlapping_gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  
  # Any better package?
  ego = enrichGO(
    gene  = gene_ETR$ENTREZID,
    keyType = "ENTREZID", 
    OrgDb   = org.Hs.eg.db,
    ont     = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  
  # barplot(ego, showCategory = 10)
  #dp = dotplot(ego, showCategory = 10)
  # cnetplot(ego, showCategory = 5)
  return(ego)
}

save_GO = function(ego, filename, dir){
  plt.list = list(dotplot(ego), cnetplot(ego, showCategory = 3))
  marrangeGrob(grobs=plt.list,ncol=1,nrow=length(plt.list)) %>% 
    ggsave(file.path(dir, filename),
           plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)
}
