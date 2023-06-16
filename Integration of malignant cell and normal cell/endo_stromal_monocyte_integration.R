#! /usr/bin/Rscript

library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)



####################################
############ Presets    
####################################

anno_file <- "adata_for_seurat_monocyte_no_AML_meta.csv"
adata_name <- "adata_for_seurat_monocyte_no_AML"
output_name <- "adata_for_seurat_monocyte_no_AML"

pc_select <- 30
res <- 0.4
batch_used <- "cancer_type"


####################################
############ Loading
####################################


Convert(paste(adata_name, ".h5ad", sep = ""), "h5seurat", overwrite = TRUE, assay = "RNA")

pbmc <- LoadH5Seurat(paste(adata_name, ".h5seurat", sep = ""), meta.data = FALSE, misc = FALSE)

anno <- read.csv(anno_file, header = TRUE, row.names = 1, check.names = F)
anno <- anno[rownames(pbmc@meta.data), ]

pbmc@meta.data$cancer_type <- anno$cancer_type
pbmc@meta.data$patient_id <- anno$study_id




####################################
############ Integration: SCTransform
####################################

# https://satijalab.org/seurat/articles/integration_introduction.html

pbmc.list <- SplitObject(pbmc, split.by = batch_used)
pbmc.list <- lapply(X = pbmc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")


####################################
############ Clustering
####################################

immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:pc_select)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:pc_select)
immune.combined <- FindClusters(immune.combined, resolution = res)

save(immune.combined, file = paste(output_name, "_sct.RData", sep = ""))



adata_name <- "adata_for_seurat_monocyte_no_AML_sct"

############ Clustering

load(paste(adata_name, ".RData", sep = ""))

# https://satijalab.org/seurat/articles/integration_introduction.html

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#immune.combined <- RunPCA(immune.combined, npcs = pc_select, verbose = FALSE)
#DimHeatmap(immune.combined, dims = 1:40, cells = 500, balanced = TRUE)
#immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:pc_select)
#immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:pc_select)

immune.combined <- FindClusters(immune.combined, resolution = 0.1)

immune.combined <- FindClusters(immune.combined, resolution = 0.2)

immune.combined <- FindClusters(immune.combined, resolution = 0.3)

immune.combined <- FindClusters(immune.combined, resolution = 0.4)

immune.combined <- FindClusters(immune.combined, resolution = 0.5)

immune.combined <- FindClusters(immune.combined, resolution = 0.6)

immune.combined <- FindClusters(immune.combined, resolution = 0.7)

immune.combined <- FindClusters(immune.combined, resolution = 0.8)

immune.combined <- FindClusters(immune.combined, resolution = 0.9)

immune.combined <- FindClusters(immune.combined, resolution = 1)

write.csv(immune.combined@meta.data, paste(adata_name, "_cluster_multiple_res.csv", sep = ""))



