library(data.table)
library(stringr)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
df = fread('/Users/jefft/Desktop/p53_project/datasets/COSMIC/V95_38_TARGETEDSCREENMUTANT.csv', na.strings = 'null')
df = as.data.frame(df)

# Identify non-silent mutations selected in different cancer types

# only consider primary tumour
df = subset(df, df$TUMOUR_ORIGIN=='primary')
# only consider mutations causing a protein change
df = subset(df, !is.na(df$HGVSP))
df = df[-grep('=', df$HGVSP), ]
# primary site must be specified
df = subset(df, df$PRIMARY_SITE != 'NS')

# only consider carcinoma ?
df = subset(df, df$PRIMARY_HISTOLOGY == 'carcinoma')
# only consider missense ?
df = subset(df, df$MUTATION_DESCRIPTION == 'Substitution - Missense')

df$PRIMARY_SITE = gsub('large_intestine', 'colon', df$PRIMARY_SITE)
df$PRIMARY_SITE = gsub('small_intestine', 'colon', df$PRIMARY_SITE)
cancer_to_analyze = c('breast', 'lung', 'stomach', 'liver', 'prostate', 
                      'ovary', 'colon', 'oesophagus', 'pancreas')
df = subset(df, df$PRIMARY_SITE %in% cancer_to_analyze)

# If only consider missense, use amino acid position to merge mutations
df['aa_pos'] = as.character(str_extract(df$MUTATION_AA,"(?<=[A-Z])[0-9]+(?=[A-Z])"))
freq_tb = table(df$aa_pos)
# Other wise, use HGVSP (per-mutation)
# freq_tb = table(df$HGVSP)

freq = as.numeric(freq_tb / nrow(df))
names(freq) = names(freq_tb)
summary(freq)
freq = sort(freq, decreasing = T)
# hist(log10(freq), breaks=50)  # set a minimum threshold of mutations?
# trh = 10**(-3)
# freq = freq[which(freq > trh)]


get_tissue_freq = function(df, cancer_site, group_by = 'aa_pos'){
  df_cancer_sub = subset(df, df$PRIMARY_SITE == cancer_site)
  freq_tb = table(df_cancer_sub[group_by])
  freq = as.numeric(freq_tb / nrow(df_cancer_sub))
  names(freq) = names(freq_tb)
  freq = sort(freq, decreasing = T)
  return(freq)
}

compare_freq = function(gt_freq, tissue_freq){
  pool = unique(c(names(gt_freq), names(tissue_freq)))
  rt = data.frame()
  for (i in pool){
    rt = rbind(rt, c(i, gt_freq[i], tissue_freq[i]))
  }
  colnames(rt) = c('aa_pos', 'Pan_cancer', 'Tissue_specific')
  rt[is.na(rt)] = 0
  rt[,2] = as.numeric(rt[,2])
  rt[,3] = as.numeric(rt[,3])
  rt['Tissue_VS_Pan'] = rt[,3] / rt[,2]
  rt = rt[order(rt[['Tissue_VS_Pan']], decreasing = T),]
  return(rt)
}

comp_pan_cancer = data.frame()
for (site in cancer_to_analyze){
  freq_tissue = get_tissue_freq(df, cancer_site = site, group_by = 'aa_pos')
  #head(freq_tissue)
  comp = compare_freq(freq, freq_tissue)
  comp_sub = subset(comp, comp$Tissue_specific >= 0.01)
  comp_sub['site'] = site
  comp_pan_cancer = rbind(comp_pan_cancer, comp_sub)
  
}
write.table(comp_pan_cancer, '/Users/jefft/Desktop/p53_project/datasets/COSMIC/tissue_p53_mut_freq.txt', 
            sep='\t', quote=F, row.names = F)
comp_pan_cancer = read.table('/Users/jefft/Desktop/p53_project/datasets/COSMIC/tissue_p53_mut_freq.txt', 
                             sep='\t', header = T)

site_to_write = 'lung'
breast_comp = subset(comp_pan_cancer, comp_pan_cancer$site==site_to_write)
breast_comp['hot-spot'] = breast_comp$aa_pos %in% names(freq)[1:6]
breast_comp['tissue_high'] = breast_comp$Tissue_VS_Pan > 1.5 & breast_comp$Tissue_specific > 0.01
hot_meta_mut = breast_comp$aa_pos[which(breast_comp$`hot-spot`)]
tissue_meta_mut = breast_comp$aa_pos[which(breast_comp$tissue_high)]
id = c(rep('hot_spot', length(hot_meta_mut)), rep('tissue_high', length(tissue_meta_mut)))
meta_muts = data.frame('meta_mut_id'=id, 'aa_pos'=c(hot_meta_mut, tissue_meta_mut))
write.table(meta_muts, '/Users/jefft/Desktop/p53_project/datasets/meta_muts/lung_COSMIC.txt',
            sep='\t', quote=F, row.names = F)

library(ggplot2)
library(plotly)
comp_pan_cancer['label'] = F
comp_pan_cancer$label[which(
  comp_pan_cancer$Tissue_VS_Pan > 2 | comp_pan_cancer$Tissue_specific > 0.2)] = T  # tissue specific hot spots
comp_pan_cancer$label[which(comp_pan_cancer$aa_pos %in% names(freq)[1:6])] = T  # hot spots

ggplot(comp_pan_cancer, 
       aes(x=Pan_cancer, y=Tissue_specific)) + theme_classic() +
  geom_point(alpha=0.5) +
  geom_abline(slope=1, intercept = 0) +
  geom_abline(slope=1.5, linetype='dotted') +
  geom_abline(slope=1/1.5, linetype='dotted') +
  geom_vline(xintercept = 0.01, linetype='dotted') +
  geom_text(data = subset(comp_pan_cancer, label), aes(label=aa_pos), size=3, nudge_x = 0.005, nudge_y = 0.005) +
  facet_wrap(~site, scales = 'free')
ggplotly()

