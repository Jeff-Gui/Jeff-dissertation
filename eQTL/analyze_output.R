library(tidyverse)
library(pheatmap)
library(patchwork)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')

gc()
output_dirs = list.files('outputs')
coll = data.frame()
for (i in output_dirs){
  if (i == 'readme.txt'){next}
  if (i == 'tcga_luad_raw_seq'){next}
  trans_eqtls = read.table(file.path('outputs', i, 'trans_eqtl_fdr005.txt'), sep='\t', header = T)
  trans_eqtls['experiment'] = i
  coll = rbind(coll, trans_eqtls)
}

coll = coll[grep('tcga', coll$experiment),]
coll = subset(coll, coll$experiment == 'tcga_brca_raw_seq')

# Positive controls
pos_control_genes = list(
  'cell_cycle' = c('CCNA2', 'CHEK1'),
  'proteosome' = c('PSMA1', 'PSMC1'),
  'nucleomtb' = c('GMPS', 'DHFR', 'IMPDH', 'AK2', 'CAD', 'DPYD',
                  'DTYMK', 'GART', 'RRM2', 'TK1', 'TYMS')
)
pos_ctr_eqtl = subset(coll, coll$gene %in% unlist(pos_control_genes) 
#                      & abs(coll$beta) > 1
                      )
## count how many has been recovered.
sm_pos = function(recovered, pos_control_genes){
  # recovered: a vector of genes
  # pos_control_genes: a list of gene groups
  ct = c()
  for (nm in names(pos_control_genes)){
    ct = c(ct, sum(recovered %in% pos_control_genes[[nm]]))
  }
  names(ct) = names(pos_control_genes)
  return(ct)
}
count_pos_ctr = pos_ctr_eqtl %>% group_by(experiment, protein_change) %>%
                        summarise(recovered_pos = sm_pos(gene, pos_control_genes))
count_pos_ctr['ctr_group'] = rep(names(pos_control_genes), nrow(count_pos_ctr) / length(pos_control_genes))
count_pos_ctr['full_pos'] = sapply(count_pos_ctr$ctr_group, function(x)return(length(pos_control_genes[[x]])))
ctr_plt = ggplot(count_pos_ctr, aes(x=ctr_group)) + theme_classic() +
  geom_bar(aes(y=full_pos), stat = 'identity', position = 'dodge', fill='gray', alpha=0.5) +
  geom_bar(aes(y=recovered_pos, fill=experiment), stat = 'identity', position = 'dodge') +
  facet_wrap(~protein_change) + labs(y='control gene count', x='category') +
  mytme + theme(strip.text = element_text(size=15), axis.text.x = element_text(angle=0))
ctr_plt
ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/eQTL', 'control_gene_plt_beta_1.pdf'),
       plot=ctr_plt, width=15,height=8,units='in',device='pdf',dpi=300)


# Discovery
mut_to_see = 'hot_spot'  # tissue_high, breast_w2016, hot_spot, p.R248Q, p.R175H
mut = subset(coll, coll$protein_change == mut_to_see)
m = matrix(0, nrow = length(unique(mut$gene)), ncol = length(unique(mut$experiment)))
m_fdr = matrix(0, nrow = length(unique(mut$gene)), ncol = length(unique(mut$experiment)))
rownames(m) = unique(mut$gene)
colnames(m) = unique(mut$experiment)
rownames(m_fdr) = unique(mut$gene)
colnames(m_fdr) = unique(mut$experiment)
for (i in unique(mut$gene)){
  sub = subset(mut, gene==i)
  for (e in 1:nrow(sub)){
    m[i,sub$experiment[e]] = abs(sub$beta[e])
    m_fdr[i,sub$experiment[e]] = abs(sub$FDR[e])
  }
}
beta = pheatmap(m, cluster_rows = F, cluster_cols = F, labels_row = '', angle_col = 45)
fdr = pheatmap(m_fdr, cluster_rows = F, cluster_cols = F, labels_row = '', angle_col = 45)


gene_nm = list()
for (i in unique(mut$experiment)){
  gene_nm[[i]] = subset(mut, mut$experiment==i)$gene
}

library(ggvenn)
library(ggplotify)
venn_plt = ggvenn(gene_nm, stroke_size = 0.5, set_name_size = 4) + 
  labs(title = mut_to_see)


overlap_up = get_intersection_eqtl(subset(mut, mut$beta > 1), group_col = 'experiment',
                                groups = c('tcga_brca_raw_seq')
                                )
overlap_down = get_intersection_eqtl(subset(mut, mut$beta < -1), group_col = 'experiment',
                                   groups = c('tcga_brca_raw_seq')
)
#overlap = get_intersection_eqtl(mut, group_col = 'experiment', 
#                                groups = c('tcga_brca_z_score')
#                                )


# Ontology
library(tidyverse)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
library(gridExtra)

do_GO = function(df){
  avg_score = aggregate(beta~gene, df, mean)
  avg_score = avg_score[order(avg_score$beta, decreasing = T),]
  overlapping_gene = avg_score$gene
  gene_ETR = bitr(overlapping_gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  #mart = useMart("ensembl", "hsapiens_gene_ensembl")
  #gene_ETR = getBM(attributes = "entrezgene_id", filters = "symbol", values = overlapping_gene, mart = mart)
  
  # Any better package?
  ego = enrichGO(
    gene  = gene_ETR$ENTREZID,
    keyType = "ENTREZID", 
    OrgDb   = org.Hs.eg.db,
    ont     = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  
  # barplot(ego, showCategory = 10)
  dp = dotplot(ego, showCategory = 10)
  # cnetplot(ego, showCategory = 5)
  return(dp)
}
up = do_GO(overlap_up)
down = do_GO(overlap_down)

if (nrow(down$data) == 0){
  plot.list = list('venn'=venn_plt, 'up_GO'=up)
} else {
  plot.list = list('venn'=venn_plt, 'up_GO'=up, 'down_GO'=down) 
}

marrangeGrob(grobs=plot.list,ncol=1,nrow=length(plot.list)) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/eQTL', paste('GO_',mut_to_see,'.pdf', sep='')),
         plot=., width=8.27,height=11.69,units='in',device='pdf',dpi=300)

# Pathway enrichment
kegg.gmt = read.gmt('/Users/jefft/Genome/c2.cp.kegg.v7.4.symbols.gmt')
geneList = avg_score$beta
names(geneList) = avg_score$gene
gsea = GSEA(geneList, TERM2GENE = kegg.gmt, pvalueCutoff = 0.05)
dotplot(gsea)
gseaplot2(gsea, 1)

# Plot representitive genes
source('load_data_cbp.R')
source('utils.R')
library(yaml)
load(file = '/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/clean_data.RData')
config_name = 'tcga_brca.yaml'

default_cfg = yaml.load_file(file.path('config', 'default.yaml'))
if (config_name != 'dafault.yaml'){
  config = yaml.load_file(file.path('config', config_name))
  config = merge_cfg(default_cfg, config)
}
eqtl_cfg = config$eQTL

eqtl_m = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
                    sample_col_name = eqtl_cfg$maf_sample_col_nm,
                    min_sample_per_snp = eqtl_cfg$min_vaf,
                    meta_mut_fp = eqtl_cfg$meta_mut)
snps = as.matrix(eqtl_m$snp)

mut = mut[order(mut$beta, decreasing = T),]
top_up = mut$gene[1:5]
snpid_to_plot = mut$snps[1]

plt.list = list()
for (top_gene in top_up){
  gene_to_plot = top_gene
  plot_df = get_single_eqtl_plot_dt(gene_to_plot, snpid_to_plot, dt, as.matrix(snps))
  plot_title = mut_to_see
  a = ggplot(plot_df, aes(x=as.logical(mut_state), y=expression)) + theme_classic() +
    geom_violin(width = 1) +
    # geom_boxplot(width = 0.05, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha=0.2) +
    labs(x=plot_title, title = strsplit(config_name, split='\\.')[[1]][1], y=gene_to_plot)
  
  if (sum(plot_df$mut_state) != sum(as.logical(plot_df$mut_state))){
    b = ggplot(subset(plot_df, mut_state != 0), aes(x=mut_state, y=expression)) + 
      theme_classic() +
      geom_point() +
      labs(x='VAF', y=gene_to_plot)
    plt.list[[paste(top_gene, '_violin')]] = a
    plt.list[[paste(top_gene, '_dot')]] = b
  } else {
    plt.list[[top_gene]] = a
  }
}

ncol = 4
nrow = ceiling(length(plt.list)/ncol)
marrangeGrob(grobs=plt.list,ncol=ncol,nrow=nrow,
             layout_matrix = matrix(1:(ncol*nrow), ncol=ncol, byrow = T)) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/eQTL', paste('Top5_',mut_to_see,'.pdf', sep='')),
         plot=., width=11.69,height=8.27,units='in',device='pdf',dpi=300)


