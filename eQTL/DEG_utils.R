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
outdir = '/Users/jefft/Desktop/p53_project/deg_experiments/TCGA_BRCA'
source('/Users/jefft/Desktop/p53_project/scripts/eQTL/utils.R')
source('/Users/jefft/Desktop/p53_project/scripts/eQTL/enrich_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
load('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData')
mut_ann = annotate_sample_mut(dt[[2]]@data)

## Laod controls ====
library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
ctrs = ctrs[-grep('wt', ctrs$Gene_annotation),]
pos_ctr = unique(ctrs$Gene)

## Annotating p53 ====
missense = mut_ann$missense
# hot_spot = colnames(snps)[which(snps['hot_spot',]!=0)]
nonsense = mut_ann$nonsense
wildtype = rownames(dt[[1]]@colData)
all_mut = unlist(mut_ann)
names(all_mut) = NULL
wildtype = wildtype[which(!wildtype %in% all_mut)]
fsv_shift = mut_ann$frameshift
expression = assay(dt[[1]][,,'RNA'])

## Run DEG ====
# degs = do_deg(nonsense, hot_spot, expression)
degs = do_deg(wildtype, all_mut, expression)
print(sum(pos_ctr %in% degs$gene))
write.table(degs, file.path(outdir, 'degs.txt'), sep='\t', row.names = F, quote = F)

degs = read.table(file.path(outdir, 'degs.txt'), header = T)
ggplot(degs, aes(x=diff, y=-log10(FDR))) +
  geom_point(size=0.5)

trh = seq(0,1,0.05)
pos_len = c()
for (i in trh){
  difgene = degs$gene[which(abs(degs$diff)>i)]
  pos_ctr_ftd = pos_ctr[which(pos_ctr %in% difgene)]
  pos_len = c(pos_len, length(pos_ctr_ftd))
}
g = ggplot(data.frame('threshold'=trh, 'ctr_gene_nm'=pos_len)) +
  geom_vline(xintercept = 0.8, color='red', size=1) +
  geom_point(aes(x=threshold, y=ctr_gene_nm)) +
  labs(x='ES threshold', y='# Control DEGs') +
  scale_x_continuous(breaks=seq(0,1,0.2)) +
  mytme
ggsave(file.path(outdir, 'filter_DEGs.pdf'),
       plot=g, width=4,height=3,units='in',device='pdf',dpi=300)
difgene = degs$gene[which(abs(degs$diff)>0.8)]

## Enrichment ====
test = do_GO(difgene, background = rownames(assay(dt[[1]][,,'RNA']))) # cell cycle markers

## Clustering ====
expr = expression[difgene, c(all_mut, wildtype)]
library(pheatmap)
ann_spl = data.frame('p53_state'=rep('wildtype', length(c(wildtype, all_mut))))
rownames(ann_spl) = c(wildtype, all_mut)
for (i in names(mut_ann)){
  ann_spl[mut_ann[[i]],'p53_state'] = i
}
table(ann_spl)
ann_spl = cbind(ann_spl, as.data.frame(dt[[1]]@colData[rownames(ann_spl),'SUBTYPE']))
colnames(ann_spl) = c('p53_state', 'subtype')
ann_spl_sub = subset(ann_spl, ann_spl$p53_state %in% c('wildtype', 'missense') & ann_spl$subtype!='')
expr_sub = expr[,rownames(ann_spl_sub)]
htmap = pheatmap(expr_sub, show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = T,
         annotation_col = ann_spl_sub, file=file.path(outdir, 'BRCA_heatmap.png'))

#extract cluster from the heatmap
clst = htmap$tree_col
out.id = cutree(clst, k=3)
ann_spl_sub$cluster = out.id
table(subset(ann_spl_sub, ann_spl_sub$p53_state != 'wildtype')$cluster)
# plot(clst)
dt[[1]]@colData$cluster = NA
dt[[1]]@colData[rownames(ann_spl_sub),'cluster'] = ann_spl_sub$cluster

# stat of missense mutation
freq_in_maf = function(cluster_nm, dt){
  maf_cst3 = rownames(ann_spl_sub)[ann_spl_sub$cluster==cluster_nm]
  maf_cst3 = subsetMaf(dt[[2]], genes = 'TP53', tsb=maf_cst3)
  table(maf_cst3@data$Protein_position)
  cst3_mut_pos = maf_cst3@data$Protein_position
  cst3_mut_freq = as.data.frame(table(cst3_mut_pos) / length(cst3_mut_pos))
  colnames(cst3_mut_freq) = c('pos', paste('frequency_', cluster_nm, sep=''))
  cst3_mut_freq[paste('count_', cluster_nm, sep='')] = table(cst3_mut_pos)
  return(cst3_mut_freq)
}
cst3_freq = freq_in_maf(cluster_nm = 3, dt)
cst2_freq = freq_in_maf(cluster_nm = 2, dt)
coll_freq = full_join(cst3_freq, cst2_freq, by='pos')
coll_freq[is.na(coll_freq)] = 0
coll_freq$sm_count = coll_freq$count_2 + coll_freq$count_3
write.table(coll_freq, file.path(outdir, 'ggplot_df.txt'), sep='\t', quote = F)
coll_freq = read.table(file.path(outdir, 'ggplot_df.txt'), sep='\t')
library(ggrepel)
coll_freq_long = gather(coll_freq %>% mutate(direction = frequency_3>frequency_2), 
                        key='gp',value='freq', c(2,4))
coll_freq_long = coll_freq_long[order(coll_freq_long$direction, coll_freq_long$freq),]
coll_freq_long$pos = as.character(coll_freq_long$pos)
coll_freq_long$pos = factor(coll_freq_long$pos, levels = unique(coll_freq_long$pos))
g = ggplot(coll_freq_long, aes(x=pos, y=freq, color=gp, group=pos)) +
  theme_classic() +
  geom_point(position = position_dodge(width = 1)) +
  geom_path(color='black', arrow = arrow(ends = 'last', length = unit(0.05, 'inches')),
            aes(linetype = direction)) +
  mytme +
  labs(x='Position', y='Frequency') +
  theme(axis.text.x = element_text(size=6, face='bold'))
g

# g = ggplot(coll_freq, aes(x=frequency_2, y=frequency_3)) +
#   geom_point(aes(size=sm_count)) +
#   geom_text_repel(aes(label=pos), color='red')
ggsave(file.path(outdir, 'freq_in_cst2_3.pdf'),
       plot=g, width=11.67,height=4,units='in',device='pdf',dpi=300)

p53_exp = as.data.frame(t(expression['TP53',]))
p53_exp['state'] = NA
p53_exp[rownames(p53_exp)[which(!rownames(p53_exp) %in% unlist(mut_ann))], 'state'] = 'wt'
p53_exp[hot_spot, 'state'] = 'exp'
p53_exp[nonsense, 'state'] = 'non'
p53_exp[fsv_shift, 'state'] = 'frameshift'
p53_exp = na.omit(p53_exp)
ggplot(p53_exp, aes(x=state, y=TP53)) +
  geom_violin()



