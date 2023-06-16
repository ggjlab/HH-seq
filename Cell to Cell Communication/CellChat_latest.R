#! /usr/bin/Rscript


#Loading

#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(NMF)
library(ggalluvial)



args <- commandArgs()

patient_id <- args[6]

n_core <- 4

#CellChat requires two user inputs: one is the gene expression data of cells, and the other is either user assigned cell labels (i.e., label-based mode) or a low-dimensional representation of the single-cell data (i.e., label-free mode). For the latter, CellChat automatically groups cells by building a shared neighbor graph based on the cell-cell distance in the low-dimensional space or the pseudotemporal trajectory space.

Convert(paste(patient_id, "_DGE_for_CellChat.h5ad", sep = ""), dest = "h5seurat", overwrite = TRUE, assay = "RNA")

pbmc <- LoadH5Seurat(paste(patient_id, "_DGE_for_CellChat.h5seurat", sep = ""), meta.data = FALSE, misc = FALSE)

#For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames. Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")

anno <- read.csv(paste(patient_id, "_DGE_for_CellChat.csv", sep = ""), header = TRUE, row.names = 1, check.names = F)

pbmc <- pbmc[, rownames(anno)]

pbmc@meta.data$ccc_anno <- anno$ccc_anno

Idents(pbmc) <- "ccc_anno"

#Part I: Data input & processing and initialization of CellChat object

cellchat <- createCellChat(object = pbmc, group.by = "ccc_anno")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

#https://github.com/sqjin/CellChat/issues/312
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = n_core) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human)

# Part II: Inference of cell-cell communication network

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#https://rstudio-pubs-static.s3.amazonaws.com/879013_4f5e5266b5824e5d8d6ddc5c357a5b37.html
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest.
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, paste(patient_id, "_net.csv", sep = ""))
write.csv(df.netp, paste(patient_id, "_netp.csv", sep = ""))

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file = paste(patient_id, "_CellChat.rds", sep = ""))











