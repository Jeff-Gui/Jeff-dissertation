library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(gridExtra)
library(enrichplot)

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

do_GO = function(df, background=NULL, ont='BP'){
  if (typeof(df)=='list'){
    avg_score = aggregate(beta~gene, df, mean)
    avg_score = avg_score[order(avg_score$beta, decreasing = T),]
    overlapping_gene = avg_score$gene
  } else {
    overlapping_gene = df
  }
  
  # gene_ETR = bitr(overlapping_gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  # if (convert_bg){
  #   BG_ETR = bitr(background, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  #   background = BG_ETR$ENTREZID
  # }
  # ego = enrichGO(gene = gene_ETR$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
  #                ont = ont, universe = background, pAdjustMethod = "BH",
  #                 pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    
  ego = enrichGO(gene = overlapping_gene, keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                  ont = ont, universe = background, pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = FALSE)
  
  # barplot(ego, showCategory = 10)
  #dp = dotplot(ego, showCategory = 10)
  # cnetplot(ego, showCategory = 5)
  return(ego)
}

save_GO = function(ego, filename, dir, cneplt=TRUE,size=NULL){
  tryCatch({
    width = 8.27
    height = 11.69
    if (cneplt){
      plt.list = list(dotplot(ego), cnetplot(ego, showCategory = 3))
      height = 11.69
    } else {
      plt.list = list(dotplot(ego))
      height = 11.69*0.4
    }
    if (!is.null(size)){
      width = size[1]
      height = size[2]
    }
    marrangeGrob(grobs=plt.list,ncol=1,nrow=length(plt.list)) %>% 
      ggsave(file.path(dir, filename),
             plot=., width=width,height=height,units='in',device='pdf',dpi=300)
  }, error = function(e){
    plt.list = list(dotplot(ego))
    marrangeGrob(grobs=plt.list,ncol=1,nrow=length(plt.list)) %>% 
      ggsave(file.path(dir, filename),
             plot=., width=width,height=height,units='in',device='pdf',dpi=300)
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
}


pcs_GO_out = function(ern, pvaluecutoff = 0.05, filename = NULL, dir = NULL, cneplt=TRUE, size=NULL){
  if (is.null(ern)){
    return(NULL)
  }
  sig_p = which(ern@result$p.adjust<0.05)
  msg = length(sig_p)
  if (msg > 0){
    if (!is.null(filename) & !is.null(dir)){save_GO(ern, filename, dir, cneplt = cneplt, size=size)}
    return(list('n_sig_term' = msg, 'result' = ern@result[sig_p,]))
  } else {
    return(list('n_sig_term' = msg, 'result' = NULL))
  }
}


load_GO_out = function(exp_home, 
                       fname = 'GO_result_no_BG_filter.RData',
                       gene_ratio_min = 0){
  load(file.path(exp_home, fname))
  go_df = data.frame()
  for (i in 1:length(go_coll)){
    nm = names(go_coll)[i]
    sub_df = go_coll[[nm]]
    sub_df[['experiment']] = nm
    rownames(sub_df) = NULL
    go_df = rbind(go_df, sub_df)
  }
  go_df['cancer'] = sapply(go_df$experiment, 
                           function(x){toupper(strsplit(x, split = '_')[[1]][2])})
  go_df['sign'] = 'neg'
  go_df[grep('pos', go_df$experiment), 'sign'] = 'pos'
  go_df['mutation'] = sapply(go_df$experiment,
                             function(x){strsplit(x, split='-')[[1]][2]})
  go_df$qvalue[which(is.na(go_df$qvalue))] = go_df$p.adjust[which(is.na(go_df$qvalue))]
  
  # QC gene ratio > 1%
  go_df['gene_ratio_val'] = sapply(go_df$GeneRatio, function(x){eval(parse(text = x))})
  
  #hist(go_df$gene_ratio_val[grep('GO', go_df$ID)], breaks=100)
  #summary(go_df$gene_ratio_val[grep('GO', go_df$ID)])
  go_df = go_df[which(go_df$gene_ratio_val > gene_ratio_min),]
  return(list('result' = result, 'df' = go_df))
}

get_DEG = function(x, gene_nm='gene', es_nm='beta', 
                   p_nm='p.adjust', p_cutoff=0.05, top_n=NULL,
                   es_cutoff=NULL){
  rt = list()
  x = subset(x, x[[p_nm]]<p_cutoff)
  if (!is.null(es_cutoff)){
    x = subset(x, abs(x[[es_nm]]) > es_cutoff)
  }
  x_pos = subset(x, x[[es_nm]]>0)
  x_neg = subset(x, x[[es_nm]]<0)
  if (is.null(top_n)){
    top_n = nrow(x)
  }
  if (nrow(x_pos)>0){
    x_pos = x_pos[order(x_pos[[es_nm]], decreasing = T),]
    x_pos = x_pos[1:min(top_n, nrow(x_pos)),]
    rownames(x_pos) = NULL
    rt[['pos']] = x_pos
  }
  if (nrow(x_neg)>0){
    x_neg = x_neg[order(abs(x_neg[[es_nm]]), decreasing = T),]
    x_neg = x_neg[1:min(top_n, nrow(x_neg)),]
    rownames(x_neg) = NULL
    rt[['neg']] = x_neg
  }
  return(rt)
}


do_GSEA = function(gene, kegg_fp = '/Users/jefft/Genome/c2.cp.kegg.v7.5.1.symbols.gmt',
                   rank_nm = 'beta', gene_nm = 'gene', pvalue=0.05){
  if (typeof(gene)=='list'){
    gene = gene[order(gene[[rank_nm]], decreasing = T),]
    gene_nm = gene[[gene_nm]]
    gene = gene[[rank_nm]]
    names(gene) = gene_nm
  }
  kegg_gmt = read.gmt(kegg_fp)
  if (gene[1] * gene[length(gene)] > 0){
    if (gene[1] > 0){
      st = 'pos'
    } else {
      st = 'neg'
    }
  } else {
    st = 'std'
  }
  gsea = GSEA(gene, TERM2GENE = kegg_gmt, pvalueCutoff = pvalue, scoreType=st,
              pAdjustMethod='fdr')
  return(gsea)
}


