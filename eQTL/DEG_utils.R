### Utils for differential gene expression analysis

do_deg = function(groupA, groupB, expression){
  # simple t test of DEG
  coll = data.frame()
  gene_nm = rownames(expression)
  for (i in 1:nrow(expression)){
    if (i%%5000==0){
      print(paste('Processed ', i, '/', nrow(expression),'.', sep=''))
    }
    a = as.numeric(expression[i, groupA])
    b = as.numeric(expression[i, groupB])
    p.value = t.test(a, b)$p.value
    dif = mean(b) - mean(a)
    coll = rbind(coll, c(gene_nm[i], dif, p.value))
  }
  colnames(coll) = c('gene', 'diff', 'p.value')
  p.adj = p.adjust(coll$p.value, method = 'fdr')
  coll[['FDR']] = p.adj
  coll = subset(coll, coll$FDR < 0.05)
  coll = coll[order(coll$diff, decreasing = T),]
  return(coll)
}

# example: compare nonsense to missense hotspot
hot_spot = colnames(snps)[which(snps['hot_spot',]!=0)]
nonsense = mut_ann$nonsense
fsv_shift = mut_ann$frameshift
expression = assay(dt[[1]][,,'RNA'])
degs = do_deg(nonsense, hot_spot, expression)

p53_exp = as.data.frame(t(expression['TP53',]))
p53_exp['state'] = NA
p53_exp[rownames(p53_exp)[which(!rownames(p53_exp) %in% unlist(mut_ann))], 'state'] = 'wt'
p53_exp[hot_spot, 'state'] = 'exp'
p53_exp[nonsense, 'state'] = 'non'
p53_exp[fsv_shift, 'state'] = 'frameshift'
p53_exp = na.omit(p53_exp)
ggplot(p53_exp, aes(x=state, y=TP53)) +
  geom_violin()



