setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult TCGA-pan_VS-wt
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'coreVScontact')
data_out = file.path(dir_home, 'data_out')
source('enrich_utils.R')
source('../overlap_utils.R')
source('../set_theme.R')

library(RcisTarget)
# featherURL = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather" 
# download.file(featherURL, destfile=basename(featherURL))

# Load data ====

load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
background = rownames(mtx)[mtx$BRCA==1]

# data motif -> TF
data(motifAnnotations_hgnc)

# input gene list

#gl = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/data_out/BRCA_migration_qtl.txt', sep='\t', header = T)
#geneLists = gl$gene[gl$peak.over.chek1==1]

gl = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/outputs/tcga_brca_raw_seq/trans_eqtl_fdr005.txt', sep='\t',header = T)
gl = gl[gl$protein_change %in% c('sandwich', 'conformation', 'contact'),]
gl = gl[order(gl$beta, decreasing = T),]
coll = list()
for (i in unique(gl$protein_change)){
  coll[[i]] = gl$gene[gl$protein_change==i & gl$beta>0]
}
um = load_controls(gen_upSet_mtx(coll))
colnames(um)
gs = rownames(um)[get_idx_condt(um, test_condition = c('contact'=T, 'peak.over.chek1'=NA,
                                     'sandwich'=NA, 'conformation'=NA, 'known.wt.up'=F, 'known.wt.down'=F))]
geneLists = gs
gl = gl[gl$gene %in% gs,]
plot(gl$beta)
geneLists = gl[1:300,'gene']

# Use database?
htft = read.table('/Users/jefft/Desktop/p53_project/datasets/TF/hTFtarget/TF-Target-information.txt', sep='\t',header = T)
tf_pool = htft[htft$target=='NFATC1',]
tf_pool = tf_pool[grep('breast',tf_pool$tissue),]
test = gl[gl$gene %in% intersect(tf_pool$TF, gs),]

# load motif based on background
dbPath = "/Users/jefft/Documents/RcisTarget/hg19-tss-centered-10kb-7species.mc9nr.feather"
# data gene -> motif
# motifRankings = importRankings(dbPath)
## data gene with background ranking; background: all genes in mutant group
motifRankings = importRankings(dbPath, columns=rownames(um))
motifRankings = reRank(motifRankings)

# Run ====
## filter geneLists out those not in motifRanking database
allgene = colnames(motifRankings@rankings)
allgene = allgene[2:length(allgene)]
geneLists = intersect(geneLists, allgene)
motifErn_gene = cisTarget(geneLists, motifRankings=motifRankings,
                          motifAnnot=motifAnnotations_hgnc)
Ern_Highconf = motifErn_gene[motifErn_gene$TF_highConf!='',]
to_rm = grep('Orthology', Ern_Highconf$TF_highConf)
if (length(to_rm)>0){
  Ern_Highconf = Ern_Highconf[-to_rm,]
}
motifErn_gene$enrichedGenes[1]

# step by step
motifs_AUC = calcAUC(geneLists, motifRankings, nCores=1)
auc = getAUC(motifs_AUC)[1,]
hist(auc,breaks=150) # should be normal





