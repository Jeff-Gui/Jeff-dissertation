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
  if (typeof(df)=='list'){
    avg_score = aggregate(beta~gene, df, mean)
    avg_score = avg_score[order(avg_score$beta, decreasing = T),]
    overlapping_gene = avg_score$gene
  } else {
    overlapping_gene = df
  }
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
  tryCatch({
    plt.list = list(dotplot(ego), cnetplot(ego, showCategory = 3))
    marrangeGrob(grobs=plt.list,ncol=1,nrow=length(plt.list)) %>% 
      ggsave(file.path(dir, filename),
             plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)
  }, error = function(e){
    plt.list = list(dotplot(ego))
    marrangeGrob(grobs=plt.list,ncol=1,nrow=length(plt.list)) %>% 
      ggsave(file.path(dir, filename),
             plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)
  })
}


load_TRUST_term2gene = function(fp = '/Users/jefft/Desktop/p53_project/datasets/TF/TRUST/trrust_rawdata.human.tsv',
                                min_gene_per_set = 20, filter_known = TRUE){
  # Load a curated set of transcription factor - gene transcription (TRUST)
  # return term2gene for enrichment analysis.
  dt = read.table(fp, sep = '\t', header = F)
  colnames(dt) = c('TF', 'target', 'sign', 'PMID')
  # hist(as.numeric(table(dt$TF)), breaks=100)
  count = table(dt$TF)
  nm_qc = names(count)[which(count > min_gene_per_set)]
  dt = subset(dt, dt$TF %in% nm_qc)
  if (filter_known) {
    dt = subset(dt, dt$sign != 'Unknown')
  }
  term = paste(dt$TF, dt$sign, sep='_')
  term2gene = data.frame('term' = term, 'gene' = dt$target)
  return(term2gene)
  
  # Testing. enrichment in BRCA eQTL genes
  # genes = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg/outputs/tcga_brca_raw_seq/trans_eqtl_fdr005.txt', 
  #                    sep='\t', header = T)
  # genes = subset(genes, genes$protein_change == 'hot_spot')
  # genes_qc = subset(genes, genes$beta > 0)
  # ern = enricher(genes_qc$gene, TERM2GENE = term2gene, 
  #                universe = hg19$geneid,
  #                pvalueCutoff = 0.05)
  # rs = ern@result
  # dotplot(ern)
}


pcs_GO_out = function(ern, pvaluecutoff = 0.05, filename = NULL, dir = NULL){
  sig_p = which(ern@result$p.adjust<0.05)
  msg = length(sig_p)
  if (msg > 0){
    if (!is.null(filename) & !is.null(dir)){save_GO(ern, filename, dir)}
    return(list('n_sig_term' = msg, 'result' = ern@result[sig_p,]))
  } else {
    return(list('n_sig_term' = msg, 'result' = NULL))
  }
}
