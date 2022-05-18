library(tidyverse)
dt_dir = '/Users/jefft/Desktop/p53_project/Science2019/RNAseq/GSE131592_RAW'
source('/Users/jefft/Desktop/Manuscript/set_theme.R')

## GSE131592
### GSM3790428...3790451: K562 DMSO <genotype> rep1~3
#### Genotype: R175H, Y220C, M237I, R248Q, R273H, R282W, KO, WT
seq_begin = 3790428
seq_end = 3790451
fls = list.files(dt_dir)
genotype = rep(c('R175H', 'Y220C', 'M237I', 'R248Q', 'R273H', 'R282W', 'KO', 'WT'),
               each = 3)
reps = rep(c('r1', 'r2', 'r3'), 8)
dt = list()
for (i in seq(seq_begin, seq_end)){
  gt = genotype[i-seq_begin+1]
  rp = reps[i-seq_begin+1]
  fn = fls[grep(i, fls)]
  mtx = read.table(file.path(dt_dir, fn), sep = '\t', skip = 4, row.names = 1)
  rnm = rownames(mtx) # assume all files have same gene order
  mtx = as.matrix(mtx[,2])  # Choose 2 as the forward stranded mapping
  dt[[paste(gt, rp, sep = '_')]] = mtx
}
dt = data.frame(dt)
rownames(dt) = rnm
# Remove 'ESNGR' genes
dt = dt[-grep('ENSGR', rownames(dt)),]
# Remove not detected
dt = dt[-which(rowSums(dt) == 0),]

# PCA
dt_scaled = t(as.matrix(scale(t(dt)))) # standardization


# Cached BioMart GRCH37.p13
mart_df = read.table('/Users/jefft/Genome/GRCh37_p13_bioMart.txt', sep='\t', header = T)
mart_df = aggregate(Transcript.length..including.UTRs.and.CDS.~Gene.stable.ID.version+Gene.name,
                    mart_df, max)
rownames(mart_df) = mart_df$Gene.stable.ID.version
symbol = mart_df[rownames(dt),2:3]
print(anyNA(symbol))
colnames(symbol) = c('gene_name', 'gene_length')
# FPKM
sample_reads = colSums(dt)
fpkm = dt
lengths = symbol$gene_length
names(lengths) = rownames(symbol)
fpkm = do.call(cbind, lapply(1:ncol(dt), function(i){
  return(fpkm[,i] / sample_reads[i] / lengths * 10^9)
}))
fpkm = data.frame(fpkm)
colnames(fpkm) = colnames(dt)

saveRDS(fpkm, file=file.path(dt_dir, 'K562_DMSO_all_fpkm.rds'))
write.table(symbol, '/Users/jefft/Desktop/p53_project/Science2019/K562_DMSO_gene_id.tsv', sep = '\t',
            quote = F, row.names = T)
saveRDS(dt, file=file.path(dt_dir, 'K562_DMSO_all.rds'))


dt = readRDS(file=file.path(dt_dir, 'K562_DMSO_all.rds'))
### edgeR
flt_deg = function(x){
  x = subset(x, x$PValue<0.05)
  x = subset(x, x$logFC>2 | x$logFC<(-2))
  return(x)
}
library(edgeR)
group = factor(rep(1:8, each=3))
y = DGEList(counts = as.matrix(dt), group = group)
keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=FALSE]
y = calcNormFactors(y, method = 'TMM') # Normalize
design = model.matrix(~group)  # Convert grouping to one-hot-like format
y = estimateDisp(y, design)
# simple pairwise comparison
et = exactTest(y, pair=c(1,2))$table
ggplot(et %>%
         mutate(up = as.numeric(logFC>2 & PValue<0.05)) %>%
         mutate(down = as.numeric(logFC<(-2) & PValue<0.05)) %>%
         mutate(state=paste(up, down, sep='-')),
       aes(x=logFC, y=-log10(PValue))) +
  geom_hline(yintercept = -log10(0.05), color='gray50', linetype='dashed') +
  geom_vline(xintercept = c(-2, 2), color='gray50', linetype='dashed') +
  geom_point(size=0.5, aes(color=state)) +
  scale_color_manual(values = c('black', 'blue', 'red')) +
  mytme
et = flt_deg(et)

# quasi-likelihood F-test with GLM approach
fit = glmQLFit(y, design)
qlf = glmQLFTest(fit)
dge_de = decideTestsDGE(qlf, adjust.method = 'BH', p.value = 0.05)
topTags(qlf)

### DESeq
library(DESeq2)
# construct the design
design = data.frame('condition'=genotype, 
                    'libType'=rep('paired-end', length(genotype)),
                    'reps'=reps)
rownames(design) = colnames(dt)
# object
dds = DESeqDataSetFromMatrix(countData = as.matrix(dt),
                             colData = design,
                             design = ~condition)
dds$condition = relevel(dds$condition, ref = 'WT')
dds = DESeq(dds)
nms = resultsNames(dds)
#res = results(dds, contrast = c('condition', 'R248Q', 'WT'))
deg = list()
for (i in 2:length(nms)){
  print(nms[i])
  deg[[i-1]] = list()
  spt = strsplit(nms[i], split = '_')[[1]]
  deg[[i-1]][['group1']] = spt[2]
  deg[[i-1]][['group2']] = spt[4]
  resLFC = as.data.frame(lfcShrink(dds, coef=nms[i]))
  resLFC = na.omit(resLFC)
  resLFC = subset(resLFC, resLFC$padj < 0.05)
  deg[[i-1]][['DEG']] = resLFC
}
deg_df = data.frame()
for (i in 1:length(deg)){
  sub_df = deg[[i]][['DEG']]
  sub_df['group1'] = deg[[i]][['group1']]
  sub_df['group2'] = deg[[i]][['group2']]
  sub_df['gene_id'] = rownames(sub_df)
  sub_df = sub_df[order(sub_df$log2FoldChange, decreasing = T),]
  rownames(sub_df) = 1:nrow(sub_df)
  deg_df = rbind(deg_df, sub_df)
}

gene_matched = mart_df[deg_df$gene_id,2:3]
deg_df['gene_name'] = gene_matched$Gene.name
write.table(deg_df, '/Users/jefft/Desktop/p53_project/Science2019/K562_DMSO_DEG.tsv',
            sep='\t', quote = F, row.names = F)

# BioMart name transfer
#library('biomaRt')
#mart = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', verbose = T)
#symbol = getBM(attributes = c('transcript_length', 'ensembl_gene_id', 'ensembl_gene_id_version'),
#               filters = 'ensembl_gene_id', mart = mart, values = gene)
#sbl = aggregate(transcript_length~ensembl_gene_id, symbol, max)


### How DEG overlap among genotypes?
library(VennDiagram)
sets = list()
for (i in unique(deg_df$group1)){
  sub_df = subset(deg_df, deg_df$group1==i)
  sets[[i]] = sub_df$gene_id
}
m = matrix(0, nrow=length(unique(deg_df$group1)), ncol=length(unique(deg_df$group1)))
rownames(m) = unique(deg_df$group1)
colnames(m) = unique(deg_df$group1)
for (i in unique(deg_df$group1)){
  rw = sets[[i]]
  rw_cp = rw
  cons = rw
  for (j in unique(deg_df$group1)){
    if (i != j){
      cl = sets[[j]]
      cons = intersect(cons, cl)
      rw_cp = rw_cp[which(!rw_cp %in% intersect(rw_cp, cl))]
      m[i, j] = length(intersect(rw, cl))
    }
  }
  m[i,i] = length(rw_cp)
}
library(circlize)
chordDiagram(mat, directional = F, transparency = 0.5)

#plotMA(res, ylim=c(-2,2))
# visualisation
# we need to transform (normalize) data (DESeq analysis uses its own, not suitable for downstream)
vsd = vst(dds, blind = FALSE)
#ntd = normTransform(dds)
# heatmap
library(pheatmap)
select = order(rowMeans(counts(dds, normalized=TRUE)),
               decreasing = T)[1:100]
df = as.data.frame(colData(dds)[,c('condition', 'reps')])
pheatmap(assay(vsd)[select,], 
         cluster_rows = F, show_colnames = T, show_rownames = F, cluster_cols = F)
# distance
sampleDists = dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists), cluster_cols = T, clustering_method = 'average')
# PCA
plotPCA(vsd, intgroup=c('condition')) + mytme


### FPKM
fpkm_log = as.matrix(log2(fpkm + 1))  # Log transfer to make normal distribution
rw_mean = rowMeans(fpkm_log)
rw_sd = sqrt(rowVars(fpkm_log))
fpkm_log_z = do.call(rbind, lapply(1:nrow(fpkm_log), function(i){
  return((fpkm_log[i,] - rw_mean[i]) / rw_sd[i])
}))
rownames(fpkm_log_z) = rownames(fpkm_log)
fpkm_log_z = na.omit(fpkm_log_z)  # Remove non-expression rows
deg_fpkm = fpkm_log_z[unique(deg_df$gene_id),]
# visualize DEGs
sig_idx = which(abs(deg_df$log2FoldChange) >= quantile(abs(deg_df$log2FoldChange),0.95))
deg_df_plot = deg_df[sig_idx,]
library(pheatmap)
down_genes = deg_df_plot$gene_id[which(deg_df_plot$log2FoldChange<=0)]
up_genes = deg_df_plot$gene_id[which(deg_df_plot$log2FoldChange>0)]
pheatmap(deg_fpkm[up_genes,], show_rownames = F, cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(colors = c("#273C99","white","#FF000C"))(100),
         border_color = NA)
genes = rownames(deg_fpkm)
names(genes) = mart_df[rownames(deg_fpkm),'Gene.name']
plot(deg_fpkm[genes['PHLDA3'],])


