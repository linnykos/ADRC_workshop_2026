# Week 3: Batch effects and integration in Seurat v5
# Kevin Z. Lin, 2026-08-09
#
# This is the code from batch-integration_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data comes from the SeuratData package. Uncomment the InstallData()
# and install_github() lines the first time you run this; they download
# a few GB into your R library and only need to be run once.

# Setting up -------------------------------------------------------------------

# install.packages("devtools")
# devtools::install_github("satijalab/seurat-data")
# devtools::install_github("satijalab/azimuth")

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(Azimuth)
library(patchwork)
library(ggplot2)

rm(list = ls())

options(future.globals.maxSize = 1e10)

# The data ---------------------------------------------------------------------

# SeuratData::InstallData("pbmcsca")

obj <- SeuratData::LoadData("pbmcsca")
obj

table(obj$Method)

head(obj@meta.data)

## A filter that is itself a batch effect
obj <- subset(obj, subset = nFeature_RNA > 1000)
obj

table(obj$Method)

## The batch effect, visible in five numbers
count_mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
dim(count_mat)
count_mat[1:5, 1:5]

tapply(obj$nCount_RNA, obj$Method, median)

# Labelling the cells with Azimuth ---------------------------------------------

obj <- Azimuth::RunAzimuth(obj, reference = "pbmcref")
obj

table(obj$predicted.celltype.l1)
length(unique(obj$predicted.celltype.l2))

## Reading the confidence, not just the label
quantile(obj$predicted.celltype.l2.score)
quantile(obj$mapping.score)

## A second opinion
table(obj$predicted.celltype.l1, obj$CellType)

# Declaring the batches --------------------------------------------------------

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj

# The Week 1 pipeline, unchanged -----------------------------------------------

obj <- Seurat::NormalizeData(obj)
obj <- Seurat::FindVariableFeatures(obj)
obj <- Seurat::ScaleData(obj)

length(Seurat::VariableFeatures(obj))
head(Seurat::VariableFeatures(obj), 20)

set.seed(10)
obj <- Seurat::RunPCA(obj, verbose = FALSE, seed.use = 10)

set.seed(10)
obj <- Seurat::RunUMAP(obj,
                       dims = 1:30,
                       reduction = "pca",
                       reduction.name = "umap.unintegrated",
                       seed.use = 10)

# Looking at the damage --------------------------------------------------------

Seurat::DimPlot(obj, reduction = "umap.unintegrated", group.by = "Method")

Seurat::DimPlot(obj, reduction = "umap.unintegrated",
                group.by = "predicted.celltype.l2")

table(obj$predicted.celltype.l1, obj$Method)

# Integrating ------------------------------------------------------------------

obj <- Seurat::IntegrateLayers(object = obj,
                               method = Seurat::CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               verbose = FALSE)
obj

set.seed(10)
obj <- Seurat::RunUMAP(obj,
                       reduction = "integrated.cca",
                       dims = 1:30,
                       reduction.name = "umap.cca",
                       seed.use = 10)

# Looking again ----------------------------------------------------------------

Seurat::DimPlot(obj, reduction = "umap.cca", group.by = "Method")

Seurat::DimPlot(obj, reduction = "umap.cca",
                group.by = "predicted.celltype.l2")

# Putting the layers back ------------------------------------------------------

obj <- SeuratObject::JoinLayers(obj)
obj

score_mat <- SeuratObject::LayerData(obj,
                                     assay = "prediction.score.celltype.l1")
dim(score_mat)
rownames(score_mat)
round(score_mat[1:5, 1:5], 3)

# The two ways to get this badly wrong -----------------------------------------

## Confounding batch with your biology
table(obj$Method, obj$Experiment)
