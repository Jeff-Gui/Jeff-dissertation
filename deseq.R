library(tidyverse)
library(limma)
setwd('/Users/jefft/Desktop/p53_project/CCLE')
dt = read.csv('RNAseq.csv', header = TRUE, na.strings = '')
p53_state = read.csv('TP53 mutations.csv', header = TRUE, na.strings='')
rownames(dt) = dt[,1]
meta = dt[,1:6]
dt = dt[,7:ncol(dt)]
dt = t(dt) #  Gene ~ Cell line, TPM log2 normalized
# dt = as.matrix(dt)  

# Gather into long form (Gene, Group (cell line), Expression)
dt = data.frame(dt)
id = as.numeric(gsub('ACH\\.','', colnames(dt)))
meta['id'] = id
write.csv(meta, 'cl_meta.csv', row.names = F)
colnames(dt) = id
dt['gene'] = rownames(dt)

# Filter lineage that contains less than 5 cell lines (including engineered)
meta = meta[order(meta$lineage_1),]
to_filter = names(table(meta$lineage_1)[which(table(meta$lineage_1)<=5)])
meta = subset(meta, ! meta$lineage_1 %in% to_filter)
dt = dt[,as.character(meta$id)]

# TPM processing, scale data?, filter low exp?
EXP_TRH = 1
mean_exp = data.frame(rowMeans(dt))
dt = dt[which(mean_exp>EXP_TRH),]
mean_exp = data.frame(mean_exp[which(mean_exp>EXP_TRH),])

### ANOVA analysis
aov_p = function(gn){
  sub_gene=as.data.frame(t(dt[gn,]))
  sub_gene['lineage'] = as.vector(meta$lineage_1)
  colnames(sub_gene)[1] = 'x'
  test = summary(aov(x~lineage, data=sub_gene))[[1]]['lineage','Pr(>F)']
  return(test)
}

p = c()
count = 0
for (i in rownames(dt)){
  if (count %% 1000 == 0){
    print(count)
  }
  count = count + 1
  p = c(p, aov_p(i))
}
p_adjust = p.adjust(p, 'BH')

mean_exp['p'] = p
mean_exp['p_adjust'] = p_adjust
invar = mean_exp[which(mean_exp$p_adjust>=0.05),]
hist(mean_exp[invar$gene,])

### Limma analysis

# Create design matrix
design_list = meta$lineage_1
design = model.matrix(~design_list)

# Fit
fit = lmFit(dt, design)
fit = eBayes(fit)
output = topTable(fit,n=Inf, p.value=0.05)
var_genes = rownames(output)
invar_genes = rownames(dt)[which(!rownames(dt) %in% var_genes)]

hist(as.numeric(dt['RD3L',]), breaks = 100)
mean_exp = data.frame(rowMeans(dt))
mean_exp['invar'] = FALSE
mean_exp[invar_genes, 'invar'] = TRUE
colnames(mean_exp)[1] = 'exp'
ggplot(mean_exp) + theme_classic() +
  geom_boxplot(aes(x=invar, y=exp))

test = mean_exp[which(mean_exp$invar & mean_exp$exp>0.05),]
boxplot(test$exp)

write.csv(mean_exp, 'mean_expression.csv', row.names = T)

