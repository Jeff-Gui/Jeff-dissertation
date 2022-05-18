library(GenomicFeatures)
library(ChIPseeker)
# Make genome annotation object
grch38_ref = makeTxDbFromGFF('/Users/jefft/Genome/GRCh38_annotation/Homo_sapiens.GRCh38.99.gff3')

# Read BED file (peaks)
peaks = readPeakFile('/Users/jefft/Desktop/p53_project/Science2019/ChIPseq/SRR9090846/tagdir/peaks.bed')

# TSS
len = 1000
promoter = getPromoters(TxDb = grch38_ref, upstream = len, downstream = len)
tagmatrix = getTagMatrix(peaks, windows = promoter)
tagHeatmap(tagmatrix, xlim = c(-len,len))

plotAvgProf(tagmatrix, xlim = c(-len, len))

tss_len = 10000
peakAnno = annotatePeak(peaks, tssRegion = c(-tss_len, 0),
                        TxDb = grch38_ref, annoDb = 'org.Hs.eg.db')
plotAnnoPie(peakAnno)

anno = peakAnno@anno@elementMetadata@listData$annotation
peak_idx = c(grep('intron 1 of', anno), grep('Promoter', anno))
anno_peak = data.frame(peakAnno@anno@elementMetadata[peak_idx, ]@listData)







