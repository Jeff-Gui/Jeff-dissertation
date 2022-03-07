library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
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

