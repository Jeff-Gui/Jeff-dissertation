library(data.table)
library(tidyverse)
library(ape)
setwd('/Users/jefft/Desktop/p53_project/CCLE/RawData')

cn = data.frame(fread('CCLE_gene_cn_21Q3.csv'))
rownames(cn) = cn$V1
cn = cn[,-1]
colnames(cn) = sapply(colnames(cn), function(x){
  return(strsplit(x, split = '\\.\\.')[[1]][1])
})

eps = data.frame(fread('CCLE_expression_21Q3.csv'))
rownames(eps) = eps$V1
eps = eps[,-1]
colnames(eps) = sapply(colnames(eps), function(x){
  return(strsplit(x, split = '\\.\\.')[[1]][1])
})

# Introduce the third variable
anu = read.csv('aneuploidy_scores_21Q3.csv', header = T, row.names = 1)
idt = read.csv('sample_info_21Q3.csv', header = T, row.names = 1, na.strings = '')
idt$lineage = gsub('engineered_', '', idt$lineage)
idt$lineage = gsub('Triple Negative Breast Cancer', 'breast', idt$lineage)
idt = idt[which(!is.na(idt$lineage)),]
idt = idt[which(!idt$lineage %in% c('unknown', 'engineered')),]
map = 1:length(unique(idt$lineage))
names(map) = unique(idt$lineage)

co_gene = intersect(colnames(eps), colnames(cn))
co_cell = intersect(intersect(rownames(eps), rownames(cn)), rownames(anu))
co_cell = intersect(intersect(rownames(eps), rownames(cn)), rownames(idt))
cn = cn[co_cell,co_gene]
eps = eps[co_cell,co_gene]
anu = anu[co_cell, 'Ploidy']
idt = map[idt[co_cell,'lineage']]
out = data.frame()
lm_out = data.frame()
for (i in co_gene){
  #sm = summary(lm(eps[,i]~cn[,i]))
  #r = sm$r.squared
  dst = as.matrix(dist(cbind(eps[,i], cn[,i])))
  dst = 1/dst
  dst[is.infinite(dst)] = 0
  #out = rbind(out, Moran.I(anu, dst))
  out = rbind(out, Moran.I(idt, dst))
  lm_model = summary(lm(cn[,i]~eps[,i]))
  mtx = c(lm_model$coefficients[,1], lm_model$r.squared)
  lm_out = rbind(lm_out, mtx)
}
out['gene'] = co_gene
colnames(lm_out) = c('intercept', 'slope', 'r.squared')
coll = cbind(out, lm_out)
exp_mean = colMeans(eps)
coll['exp_mean'] = exp_mean
coll_sub = subset(coll, coll$exp_mean>1)
ggplot(coll_sub, aes(x=observed, y=r.squared)) + theme_classic() +
  geom_point(aes(color=exp_mean)) +
  geom_vline(xintercept = 0) +
  labs(x='Moran\'s I')

dt = data.frame('Gene' = colnames(eps), 'Rsq' = rs)
dt = dt[order(dt$Rsq, decreasing = T),]

plt.cor = function(x, z){
  df = data.frame('x'=eps[,x], 'y'=cn[,x], 'z'=z)
  ggplot(df, aes(x=x,y=y,color=z))+theme_classic() +
    geom_point(aes(alpha=0.5)) +
    labs(x='Expression', y='Copy Number') + 
    #scale_color_gradient2(high = 'red', low = 'blue')
    scale_x_continuous(expand = c(0,0)) +
    facet_wrap(~z, ncol = 10)
}
plt.cor('TSPAN6', as.factor(idt))


