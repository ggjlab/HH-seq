#! /usr/bin/Rscript

library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

pbmc <- LoadH5Seurat("all_obs_no_glioma_raw.h5seurat", meta.data = FALSE, misc = FALSE)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")

anno <- read.csv("all_obs_anno.csv", header = TRUE, row.names = 1, check.names = F)
anno <- anno[rownames(pbmc@meta.data), ]

pbmc@meta.data$cancer_type <- anno$cancer_type
pbmc@meta.data$if_cancer <- anno$if_cancer

Idents(pbmc) <- "cancer_type"

BRCA.markers <- FindConservedMarkers(pbmc, ident.1 = 'BRCA', grouping.var = "if_cancer")
COAD_READ.markers <- FindConservedMarkers(pbmc, ident.1 = 'COAD_READ', grouping.var = "if_cancer")
ESCA.markers <- FindConservedMarkers(pbmc, ident.1 = 'ESCA', grouping.var = "if_cancer")
HCC.markers <- FindConservedMarkers(pbmc, ident.1 = 'HCC', grouping.var = "if_cancer")
ICC.markers <- FindConservedMarkers(pbmc, ident.1 = 'ICC', grouping.var = "if_cancer")
LUAD.markers <- FindConservedMarkers(pbmc, ident.1 = 'LUAD', grouping.var = "if_cancer")
STAD.markers <- FindConservedMarkers(pbmc, ident.1 = 'STAD', grouping.var = "if_cancer")

write.csv(BRCA.markers, "BRCA.markers_no_brain.csv")
write.csv(COAD_READ.markers, "COAD_READ.markers_no_brain.csv")
write.csv(ESCA.markers, "ESCA.markers_no_brain.csv")
write.csv(HCC.markers, "HCC.markers_no_brain.csv")
write.csv(ICC.markers, "ICC.markers_no_brain.csv")
write.csv(LUAD.markers, "LUAD.markers_no_brain.csv")
write.csv(STAD.markers, "STAD.markers_no_brain.csv")


