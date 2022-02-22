library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')

gene_ann = select(org.Hs.eg.db,
                  columns = c('SYMBOL', 'GO'),
                  keys = keys(org.Hs.eg.db, keytype = 'ENTREZID'))

metabric_res = read.table('outputs/metabric/trans_eqtl_fdr005.txt', header = T)
gene_ann_ftd = subset(gene_ann, gene_ann$SYMBOL %in% metabric_res$gene)
gene_ann_ftd = subset(gene_ann_ftd, gene_ann_ftd$ONTOLOGY == 'BP') 
# BP: biological process; MF: molecular function; CC: cellular compartment

# gene_ann_ftd = subset(gene_ann_ftd, gene_ann_ftd$EVIDENCE == 'IDA') # direct process
gene_ann_ftd['GO_term'] = NA
for (i in 1:nrow(gene_ann_ftd)){
  gene_ann_ftd$GO_term[i] = GOTERM[[gene_ann_ftd$GO[i]]]@Term
}

