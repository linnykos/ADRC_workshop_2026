# https://satijalab.org/seurat/articles/pbmc3k_tutorial

rm(list=ls())
library(Seurat)
library(scCustomize)
library(EnhancedVolcano)

pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

set.seed(10)
# Let's do a simple preprocessing of pbmc
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, 
                                   project = "pbmc3k")
pbmc <- Seurat::NormalizeData(pbmc)
pbmc <- Seurat::FindVariableFeatures(pbmc, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
pbmc <- Seurat::ScaleData(pbmc)
pbmc <- Seurat::RunPCA(pbmc, 
                       features = Seurat::VariableFeatures(pbmc),
                       verbose = FALSE)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)

############

# Now for clustering
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
colnames(pbmc@meta.data)

scCustomize::DimPlot_scCustom(pbmc,
                              reduction = "umap",
                              group.by = "RNA_snn_res.0.5")

####

pbmc <- Seurat::FindClusters(pbmc, resolution = 0.1)
scCustomize::DimPlot_scCustom(pbmc,
                              reduction = "umap",
                              group.by = "RNA_snn_res.0.1")

########################

# https://satijalab.org/seurat/articles/de_vignette
Seurat::Idents(pbmc) <- "RNA_snn_res.0.1"
de_res <- Seurat::FindMarkers(pbmc,
                              ident.1 = "0", 
                              ident.2 = "1",
                              test.use = "wilcox")
head(de_res)

EnhancedVolcano::EnhancedVolcano(de_res,
                                 lab = rownames(de_res),
                                 x = "avg_log2FC",
                                 y = "p_val")

##################

# Let's find marker genes
Seurat::Idents(pbmc) <- "RNA_snn_res.0.1"
marker_res <- Seurat::FindMarkers(pbmc,
                                  ident.1 = "0", 
                                  ident.2 = NULL,
                                  test.use = "wilcox",
                                  only.pos = TRUE)
marker_res <- marker_res[marker_res$p_val_adj <= 1e-5,]
marker_res <- marker_res[order(marker_res$avg_log2FC, decreasing = TRUE),]
head(marker_res)

Seurat::DotPlot(pbmc,
                features = rownames(marker_res)[1:6],
                group.by = "RNA_snn_res.0.1")

Seurat::DotPlot(pbmc,
                features = rownames(marker_res)[51:56],
                group.by = "RNA_snn_res.0.1")

