# Ref: http://www.bio-info-trainee.com/2535.html
library(WGCNA)
library(Matrix)
library(yaml)
setwd('/Users/jefft/Desktop/p53_project')
source('scripts/eQTL/utils.R')
source('scripts/eQTL/load_data_cbp.R')

config_name = 'tcga_brca.yaml'
dt_name = 'BRCA-TCGA'

default_cfg_name = 'default.yaml'
TOM_name = 'WGCNA-TOM'
file_out_name = 'WGCNA_network.RData'
export_cyt_w_trh = 0.1
MM_trh = 0.8
GS_trh = 0.05
saveTOM = TRUE
use_cache = TRUE
run = FALSE
if (use_cache){
  saveTOM = FALSE
}
pick_soft_trh = FALSE

default_cfg = yaml.load_file(file.path('scripts/WGCNA/config', default_cfg_name))
if (config_name != default_cfg_name){
  config = yaml.load_file(file.path('scripts/WGCNA/config', config_name))
  config = merge_cfg(default_cfg, config)
}
config = pcs_cfg(config)
dt_cfg = config$dataset
eqtl_cfg = config$eQTL
preprocess_cfg = config$preprocess
dt_home = dirname(dt_cfg$dataset_home)

if (cache){
  load(file.path(dt_home, 'clean_data_WGCNA.RData'))
} else {
  refresh_log = TRUE
  ## Handle log
  logpath = file.path(config$output, 'log.txt')
  if (file.exists(logpath) & refresh_log){
    file.remove(logpath)
  }
  basicConfig(level = 'FINEST')
  addHandler(writeToFile, file=logpath, level='DEBUG')
  loginfo('Loading data...', logger = 'main')
  
  dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                     exp_nm = dt_cfg$exp_nm, 
                     cna_nm = dt_cfg$cna_nm, 
                     mut_nm = dt_cfg$mut_nm,
                     sample_meta_nm = dt_cfg$sample_meta_nm,
                     patient_meta_nm = dt_cfg$patient_meta_nm,
                     case_complete_nm = dt_cfg$case_complete_nm,
                     case_list_dir_nm = dt_cfg$case_list_dir_nm,
                     na.str = dt_cfg$na.str,
                     gene_col_nm = dt_cfg$gene_col_nm,
                     normalize_genes = preprocess_cfg$norm_gene,
                     quantile_norm = preprocess_cfg$quantile,
                     z_score = preprocess_cfg$z_score,
                     real_z_score = preprocess_cfg$z_score_real,
                     has_log_ed = dt_cfg$has_log_ed,
                     rm_low_expr_gene = as.numeric(preprocess_cfg$rm_low_expr_gene),
                     diag_out = config$output,
                     filter_protein_coding = eqtl_cfg$genepos,
                     diploid_norm = preprocess_cfg$diploid_norm)
  gc()
  save(dt, file = file.path(dirname(dt_cfg$dataset_home), 'clean_data_WGCNA.RData'))
}

p53_ann = annotate_sample_mut(dt[[2]]@data)
dt[[1]]@colData[['p53_state']] = 'wildtype'
for (i in names(p53_ann)){
  dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
}
b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, mode = 'position-tcga', 
        samples = rownames(dt[[1]]@colData), snp_list = list(c(273, 248)))
dt[[1]]@colData[rownames(b_m)[which(b_m==1)], 'p53_state'] = 'MS-contact'
b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, mode = 'position-tcga', 
                                samples = rownames(dt[[1]]@colData), snp_list = list(c(175,245,249,282)))
dt[[1]]@colData[rownames(b_m)[which(b_m==1)], 'p53_state'] = 'MS-conform'
dt[[1]]@colData[['genotype']] = dt[[1]]@colData$p53_state


if (run){
  dataExpr = assay(dt[[1]][,,'RNA'])
  allowWGCNAThreads()  # parallelize processing
  colData = as.data.frame(dt[[1]]@colData)
  type = 'signed'
  corType = 'pearson'
  corFnc = ifelse(corType == 'pearson', cor, bicor)
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  
  # do not do Median Absolute Deviation because it is on Rank-normalized data
  # m.mad = apply(dataExpr,1,mad)
  # dataExprVar = dataExpr[which(m.mad > 
  #                             max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
  # dataExprVar = as.matrix(dataExprVar)
  
  dataExpr = as.data.frame(t(dataExpr))
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  dim(dataExpr)
  
  if (pick_soft_trh){
    ###======= Pick parameters ========
    ## Soft-threshold: connectivity threshold of the nodes, for better singed network
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                            networkType=type, verbose=5)
    # par(mfrow = c(1,2))
    # cex1 = 0.9
    # 
    # # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
    # # 网络越符合无标度特征 (non-scale) -> better
    # plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    #      xlab="Soft Threshold (power)",
    #      ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    #      main = paste("Scale independence"))
    # text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    #      labels=powers,cex=cex1,col="red")
    # # 筛选标准。R-square=0.85
    # abline(h=0.85,col="red")
    # 
    # # Soft threshold与平均连通性
    # plot(sft$fitIndices[,1], sft$fitIndices[,5],
    #      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    #      main = paste("Mean connectivity"))
    # text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
    #      cex=cex1, col="red")
    power = sft$powerEstimate
    print(power)
  } else {
    if (type == 'signed'){
      power = 12  # or 12?
    } else {
      power = 7
    }
  }
  
  
  ###=========== Construct the network =============
  cor = WGCNA::cor
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=saveTOM, corType = corType, 
                         maxPOutliers=maxPOutliers,
                         saveTOMFileBase = file.path(config$output, TOM_name),
                         verbose = 3)
  cor = stats::cor
  table(net$colors)
  mergedColors = labels2colors(net$colors)
  # plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  #                     "Module colors",
  #                     dendroLabels = FALSE, hang = 0.03,
  #                     addGuide = TRUE, guideHang = 0.05)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  # MEs = net$MEs;
  geneTree = net$dendrograms[[1]]
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  # Model membership (MM)
  datKME=signedKME(dataExpr, MEs, outputColumnName="MM.")
  gene_color = list()
  probes = colnames(dataExpr)
  for (i in unique(moduleColors)){
    mmn = paste('MM.',i, sep='')
    gene_color[[i]] = probes[which(moduleColors==i & datKME[,mmn]>MM_trh)]
  }
  save(MEs, datKME, moduleLabels, moduleColors, geneTree, net, gene_color,
       file = file.path(config$output, file_out_name))
  gc()
} else {
  load(file.path(config$output, file_out_name))
}

# VIS network
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

# Module enrichment in clusters
colData = as.data.frame(dt[[1]]@colData)
colData$genotype = as.factor(colData$genotype)

# stratify by mutants
em = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
           sample_col_name = eqtl_cfg$maf_sample_col_nm,
           min_sample_per_snp = eqtl_cfg$min_vaf,
           meta_mut_fp = eqtl_cfg$meta_mut,
           mode = eqtl_cfg$mode)
lu = em[[3]]
em = as.matrix(em[[2]])
pc = sapply(rownames(em), function(x){return(get_var_info_from_maf(dt[[2]]@data, x, lu))})
pc = unlist(pc)
if (length(pc)>0){
  for (i in 1:length(pc)){
    rownames(em)[which(rownames(em) == names(pc)[i])] = pc[i]
  }
}
design = as.data.frame(t(em))
design['wildtype'] = 0
design[rownames(colData)[which(colData$genotype=='wildtype')],
       'wildtype']=1
design['nonsense'] = 0
design[rownames(colData)[which(colData$genotype=='nonsense')],
       'nonsense'] = 1

# stratify by p53 genotype
# design = model.matrix(~0 + colData$genotype)
# colnames(design) = levels(colData$genotype)
# rownames(design) = rownames(colData)

nSamples = dim(design)[1]
design = as.matrix(design)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
pdf(file = file.path(config$output, 'heatmap.pdf'), width=8, height=6)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design),
               yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE,
               colors = greenWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = T, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

## GO module genes and control gene distribution
source('scripts/eQTL/enrich_utils.R')
library(readxl) # load control genes
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = unique(ctrs$Gene)

coll_rs = data.frame()
coll_GO = list()
probes = as.character(rownames(assay(dt[[1]][,,'RNA'])))
for (nm in names(gene_color)){
  gene_in_color = probes[which(moduleColors==nm)]
  ego = do_GO(gene_color[[nm]])
  ego_rs = ego@result
  po = pcs_GO_out(ego, pvaluecutoff = 0.05, 
                  filename = paste(dt_name, '_' ,nm, '.pdf', sep=''), 
                  dir = file.path(config$output, 'GO'))
  coll_rs = rbind(coll_rs, c(length(gene_in_color), length(gene_color[[nm]]), 
                             po$n_sig_term, sum(ctrs %in% gene_in_color)))
  coll_GO[[nm]] = po$result
}
colnames(coll_rs) = c('Size', 'Size_filtered', 'GO_term', 'Control_gene')
coll_rs[['color']] = names(gene_color)
save(coll_GO, file = file.path(config$output, 'GO_result.RData'))
write.table(coll_rs, file.path(config$output, 'GO_ctr_summary.txt'), sep='\t', quote=F, row.names=F)

# Check overlap with control and genes discovered in eQTL
tb = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg/outputs/tcga_brca_raw_seq/trans_eqtl_fdr005.txt', sep='\t', header = T)
ctrs[which(ctrs %in% probes[which(moduleColors=='brown')])]
100 * table(moduleColors[which(probes %in% tb$gene)]) / length(unique(tb$gene))
table(moduleColors[which(probes %in% c('MYBL2', 'DPH2', 'DEPDC7', 'LFNG', 'BTG2'))])
rownames(datKME)[order(datKME$MM.purple, decreasing = T)][1:10]

## Identify hub genes
### Gene significance
trait_of_interest = 'hot_spot_contact'
model_of_interest = 'grey'
mmn = paste('MM.',model_of_interest, sep='')

design = as.data.frame(design)
hot_spot = as.data.frame(as.numeric(design[[trait_of_interest]]>0))
names(hot_spot) = trait_of_interest
geneTraitSignificance = as.data.frame(cor(t(assay(dt[[1]][,,'RNA'])), hot_spot, use = "p"))
names(geneTraitSignificance) = paste("GS.", names(hot_spot), sep="")
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(GSPvalue) = paste("p.GS.", names(hot_spot), sep="")
GS_MM = cbind(geneTraitSignificance, datKME[rownames(geneTraitSignificance), mmn])
colnames(GS_MM) = c('GS', 'MM')
ggplot(GS_MM) + theme_classic() +
  geom_point(aes(x=MM, y=GS), size=0.5)

filter = abs(GS_MM$MM) > GS_trh & abs(GS_MM$GS) > MM_trh
genes = rownames(GS_MM)[filter]
beta = GS_MM$MM[filter]
# ego = do_GO(data.frame('gene'=genes, 'beta'=beta))
# dotplot(ego)
# cnetplot(ego)

#======================================================

## Visualise network to cytoscape
probes = rownames(assay(dt[[1]][,,'RNA']))
# Recalculate topological overlap if not cached. VERY SLOW!!!
if (use_cache){
  load(file.path(config$output, paste(TOM_name, 'block.1.RData', sep='-')))
  # TOM = as.matrix(TOM)
} else {
  TOM = TOMsimilarityFromExpr(t(assay(dt[[1]][,,'RNA'])), power = power)
  save(TOM, file = file.path(config$output, paste(TOM_name, 'block.1.RData', sep='-')))
  gc()
}
# Select module color
coll_sign = list()
for (md in c('pink', 'brown', 'grey')){
  module = md
  # Select module probes
  inModule = (moduleColors==module)
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = export_cyt_w_trh,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  )
  ced = cyt$edgeData %>% mutate('color'=md)
  rownames(ced) = 1:nrow(ced)
  coll_sign[[md]] = ced
}
gc()
coll_sign = Reduce(rbind, coll_sign)
coll_sign = coll_sign[,which(!colnames(coll_sign) %in% c('fromAltName','toAltName'))]
write.table(coll_sign, paste('WGCNA/outputs/', dt_name, '.csv', sep=''), quote = F, row.names = F, sep=',')
remove(TOM)
