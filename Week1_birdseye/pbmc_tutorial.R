# Week 1: A bird's-eye view of a single-cell workflow
# Kevin Z. Lin, 2026-07-27
#
# This is the code from pbmc_tutorial.rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# EDIT THIS PATH to point at your own filtered_feature_bc_matrix folder:
data_dir <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix"

# Setting up the Seurat object ------------------------------------------------

rm(list = setdiff(ls(), "data_dir"))
library(Seurat)
pbmc.data <- Seurat::Read10X(data.dir = data_dir)

pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data,
                                   project = "pbmc3k")
pbmc

head(pbmc@meta.data)

count_mat <- SeuratObject::LayerData(pbmc,
                                     layer = "counts",
                                     assay = "RNA")
count_mat[42:44,14:16]
head(colSums(count_mat))
head(pbmc$nCount_RNA)

SeuratObject::Features(pbmc)[1:5]

Seurat::Cells(pbmc)[1:5]

# Quality control -------------------------------------------------------------

## Mitochondrial fraction
gene_vec <- SeuratObject::Features(pbmc)
gene_vec[grep("^MT-", gene_vec)]

pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc,
                                                     pattern = "^MT-")

## Ribosomal fraction
gene_vec[grep("^RPS", gene_vec)]
pbmc[["percent.rb"]] <- Seurat::PercentageFeatureSet(pbmc,
                                                     pattern = "^RPS")

## Cell cycle scores
pbmc <- Seurat::NormalizeData(pbmc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- Seurat::CellCycleScoring(pbmc,
                                 s.features = s.genes,
                                 g2m.features = g2m.genes)

head(pbmc@meta.data)

## Visualizing the QC metrics
Seurat::VlnPlot(pbmc, features = c("nFeature_RNA",  "nCount_RNA", "percent.mt", "percent.rb"),
                ncol = 4)

## Filtering
pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization ---------------------------------------------------------------

pbmc <- Seurat::NormalizeData(pbmc)

# Feature selection -----------------------------------------------------------

pbmc <- Seurat::FindVariableFeatures(pbmc,
                                     selection.method = "vst",
                                     nfeatures = 2000)

top10 <- head(Seurat::VariableFeatures(pbmc), 10)
plot1 <- Seurat::VariableFeaturePlot(pbmc)
plot1 <- Seurat::LabelPoints(plot = plot1,
                             points = top10,
                             repel = TRUE)
plot1

# Dimension reduction ---------------------------------------------------------

## Scaling
pbmc <- Seurat::ScaleData(pbmc)

## PCA
pbmc <- Seurat::RunPCA(pbmc,
                       features = Seurat::VariableFeatures(object = pbmc),
                       verbose = FALSE)

pbmc

# Clustering ------------------------------------------------------------------

pbmc <- Seurat::FindNeighbors(pbmc,
                              dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc,
                             resolution = 0.5)

head(pbmc@meta.data)

# Visualization ---------------------------------------------------------------

pbmc <- Seurat::RunUMAP(pbmc,
                        dims = 1:10)
plot1 <- Seurat::DimPlot(pbmc,
                reduction = "umap",
                group.by = "RNA_snn_res.0.5")
plot1

pbmc
