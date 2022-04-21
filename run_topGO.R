# sampleGOdata = new("topGOdata",
#                     description = "Simple session",
#                     ontology = "BP",
#                     allGenes = geneList[1:50],
#                     nodeSize = 10, ## 删去少于10个注释基因的GO term
#                     annot = annFUN.db) ## maps genes identifiers to GO terms from the affyLib object)

# extract term2gene from database
library(org.Hs.eg.db)
library(tidyverse)
test = select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = 'GO'),
       columns = c('SYMBOL', 'ONTOLOGY'), keytype = 'GO')
# Can filter evidence level
print(unique(test$EVIDENCE))
test = test[test$EVIDENCE %in% c(), ]

dummy = function(x){return(paste(x,collapse = '~'))}
go2term = test %>% group_by(GO) %>% summarise(genes=dummy(SYMBOL))
go2term_lst = list()
for (i in 1:nrow(go2term)){
  go2term_lst[[go2term$GO[i]]] = strsplit(go2term$genes[i], split='~')[[1]]
}
gene2go = inverseList(go2term_lst)
# check gene symbol is consistent
length(intersect(unique(test$SYMBOL), names(gene2go)))
length(unique(test$SYMBOL))
# archive list
save(gene2go, file='/Users/jefft/Desktop/p53_project/datasets/Gene2GO_HsDbv3.14.0.RData')

# Run topGO ====
library(topGO)
load('/Users/jefft/Desktop/p53_project/datasets/Gene2GO_HsDbv3.14.0.RData')
df = read.table('/Users/jefft/Desktop/p53_project/scripts/eQTL/hg38_gene_table_autosome.tsv', header = T)
bg = df$geneid # background
GOI = c()      # test set
geneList = rep(0,length(bg))
names(geneList) = bg
geneList[GOI] = 1
geneList = factor(geneList, levels=c(0,1))
GOdata = new("topGOdata", ontology = "BP", allGenes = geneList,
             annot = annFUN.gene2GO, gene2GO = gene2go)

resultFisher = runTest(GOdata, algorithm = "elim", statistic = "fisher")
allRes = GenTable(GOdata,classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 10)

