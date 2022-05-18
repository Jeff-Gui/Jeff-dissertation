library(tidyverse)
library(data.table)
setwd('/Users/jefft/Desktop/p53_project/CCLE')
meta = as.data.frame(fread('RawData/sample_info_21Q3.csv', na.strings = ''))
mut = as.data.frame(fread('RawData/CCLE_mutations_21Q3.csv', na.strings = ''))
mut['mut_uid'] = paste(mut$Reference_Allele, mut$Start_position, mut$Tumor_Seq_Allele1, sep='-')

tpm = as.data.frame(fread('RawData/CCLE_expression_21Q3.csv', na.strings=''))
# Filter with sample having expression data
meta = meta[which(meta$DepMap_ID %in% tpm$V1),] # read tpm from expression
mut = mut[which(mut$DepMap_ID %in% tpm$V1),]

# Quality control of variants suitable for eQTL
gene_qc = function(mut, genes=NULL){
  if (!is.null(genes)){
    mut = subset(mut, mut$Hugo_Symbol %in% genes)
  }
  SAMPLE_CUTOFF = 10
  SNP_ONLY = TRUE
  MISSENSE_ONLY = TRUE
  NON_SILENT = FALSE
  if (SNP_ONLY){
    mut = subset(mut, mut$Variant_Type=='SNP')
  }
  if (MISSENSE_ONLY){
    mut = subset(mut, mut$Variant_Classification == 'Missense_Mutation')
  }
  if (NON_SILENT){
    mut = subset(mut, !is.na(mut$Protein_Change))
  }
  # TODO MAF 1% cutoff
  all_snp = mut$mut_uid
  ct_snp = table(all_snp)
  ct_snp = ct_snp[which(ct_snp>=SAMPLE_CUTOFF)]
  mut = subset(mut, mut$mut_uid %in% names(ct_snp))
  mut = mut[order(mut$Hugo_Symbol, mut$Start_position),]
  return(mut)
}

# Generate the binary state of gene mutation for each cell line.
gene_mut_bin = function(meta, mut, genes=NULL){
  # Assume samples in meta but not mutation is wt full
  if (is.null(genes)){
    genes = unique(mut$Hugo_Symbol)
  }
  mut = subset(mut, mut$Hugo_Symbol %in% genes)
  m = matrix(data=0, nrow=length(unique(meta$DepMap_ID)), ncol=length(genes))
  rownames(m) = unique(meta$DepMap_ID)
  colnames(m) = genes
  for (g in genes){
    m[unique(subset(mut, mut$Hugo_Symbol==g)$DepMap_ID),g] = 1
  }
  return(m)
}

genes = c('TP53')
out = gene_qc(mut, genes=genes)
write.table(out, 'frequent_variant.tsv', sep='\t', quote = FALSE)
flt_meta = subset(meta, meta$DepMap_ID %in% out$DepMap_ID)
mut_bin = gene_mut_bin(meta, mut, genes=unique(out$Hugo_Symbol))
write.table(mut_bin, 'frequet_variant_bin.tsv', sep='\t', quote = FALSE)

var_count = aggregate(cDNA_Change~Hugo_Symbol, data=out, FUN=function(x){
  length(unique(x))}
)
var_count = var_count[order(var_count$cDNA_Change, decreasing = T),]
head(var_count)



