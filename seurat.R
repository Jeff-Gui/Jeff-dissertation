library(Seurat)
library(tidyverse)

dt = read.csv('/Users/jefft/Desktop/p53_project/download.csv', header = TRUE, na.strings = '')
meta = dt[,1:6]
rownames(meta) = meta[,1]
meta = meta[,2:ncol(meta)]

dt = dt[,7:ncol(dt)]
rownames(dt) = meta[,1]
dim(dt)
dt = t(dt)
dt = as.matrix(dt)
dt = CreateSeuratObject(dt)
dt@meta.data = cbind(dt@meta.data, meta)
gc()

dt = ScaleData(dt)
dt = FindVariableFeatures(dt, selection.method = "vst", nfeatures = 1500)
fts = VariableFeatures(dt)

# PCA
dt = RunPCA(dt, features = fts)
ElbowPlot(dt)

# Clustering
dt = FindNeighbors(dt, dims = 1:15)
dt = FindClusters(dt, resolution = 2)

dt = RunUMAP(dt, dims = 1:15)
DimPlot(dt)
DimPlot(dt, group.by = 'lineage_1')

