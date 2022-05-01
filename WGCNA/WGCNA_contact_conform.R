# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
library(WGCNA)
library(tidyverse)
library(RColorBrewer)
setwd('/Users/jefft/Desktop/p53_project')
source('scripts/eQTL/utils.R')
source('scripts/set_theme.R')

# load dataset ====
dt = load_clean_data('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data_noRankNorm.RData',
                ann_bin_mut_list = c('contact','conformation','sandwich'))
# load WGCNA result
load('scripts/WGCNA/outputs/tcga_brca_raw_seq/WGCNA_network.RData')
coldt = dt[[1]]@colData
ME_p53 = cbind(MEs, coldt[rownames(MEs), c('p53_state', 'has_contact')])
ggplot(ME_p53) +
  geom_point(aes(x=MEblack, y=MEpurple, color=p53_state)) +
  facet_wrap(~p53_state)
ggplot(ME_p53) +
  geom_boxplot(aes(x=p53_state, y=MEpurple))

# design[design>0] = 1 binary will reduce the power!
# get VAF information
mean_vafs = subset(dt[[2]]@data, Hugo_Symbol=='TP53') %>% 
  group_by(Tumor_Sample_Barcode) %>% summarise(mean_VAF = mean(VAF))
mean_vafs = as.data.frame(mean_vafs)
rownames(mean_vafs) = mean_vafs$Tumor_Sample_Barcode

design = get_binary_SNP_m_from_maf(dt[[2]]@data, 
              snp_list = list(c(248),c(273),c(248,273),
                              c(175),c(249),c(282),c(245),c(175,249,282,245),
                              c(245,248,249,282,273,175)),
              samples = rownames(MEs), mode = 'position',
              protein_change_col = 'HGVSp_Short', code_with_VAF = T)
too_few_col = which(colSums(design>0)<5)
if (length(too_few_col) > 0){
  design = design[,-too_few_col]
}
# include wiltype, missense, nonsense traits
design = as.data.frame(design)
for (i in c('nonsense', 'missense')){
  spl = rownames(coldt)[which(coldt$p53_state==i)]
  design[spl,i] = mean_vafs[spl, 'mean_VAF']
}
design[is.na(design)] = 0
design['wildtype'] = 1
design[coldt$p53_state != 'wildtype', 'wildtype'] = 0

moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitCor_long = gather(as.data.frame(moduleTraitCor) %>% mutate(ME=rownames(moduleTraitCor)), 
                             key='group', value='cor', 1:ncol(moduleTraitCor))
moduleTraitP_long = gather(as.data.frame(moduleTraitPvalue) %>% mutate(ME=rownames(moduleTraitPvalue)), 
                             key='group', value='pvalue', 1:ncol(moduleTraitPvalue))
moduleTraitCor_long['pvalue'] = moduleTraitP_long$pvalue
# df_plt = moduleTraitCor_long[-which(moduleTraitCor_long$pvalue>0.05),]
df_plt = moduleTraitCor_long
df_plt$ME = gsub('ME','', df_plt$ME)
df_plt$plabel = ''
for (i in 1:nrow(df_plt)){
  if (df_plt$pvalue[i] < 0.05){
    df_plt$plabel[i] = '*' 
  }
  if (df_plt$pvalue[i] < 0.01){
    df_plt$plabel[i] = '**'
  }
  if (df_plt$pvalue[i] < 0.001){
    df_plt$plabel[i] = '***'
  }
}
x_labels = c('R273', 'R248', 'R273, R248', 'R175', 'G245', 'R282',
             'R175, G245, R249, R282', 'All MS hotspots', 'Other MS', 'Nonsense', 'Wildtype')
names(x_labels) = c('273','248','248_273','175','245','282',
                    '175_249_282_245', '245_248_249_282_273_175', 'missense','nonsense', 'wildtype')
df_plt$group = x_labels[df_plt$group]
df_plt$group = factor(df_plt$group, levels=x_labels)
myPalette = colorRampPalette(c("royalblue","white", "coral"))
sc = scale_fill_gradient2(low = 'royalblue', mid = 'white', high = 'coral')
ggplot(df_plt, aes(x=group,y=ME)) + mytme +
  geom_tile(aes(fill=cor), color='grey') +
  geom_text(aes(label=plabel), color='black', vjust=1) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1),
        panel.border = element_rect(color='black', size=1, fill=NA),
        panel.grid.major = element_line(color='grey',linetype='dotted')) +
  sc


#================ original heatmap ===================
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
# pdf(file = file.path(config$output, 'heatmap.pdf'), width=8, height=6)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design),
               yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = T, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))


