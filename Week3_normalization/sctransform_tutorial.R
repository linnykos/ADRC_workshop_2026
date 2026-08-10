# Week 3: Normalizing counts and selecting genes with SCTransform
# Kevin Z. Lin, 2026-08-09
#
# This is the code from sctransform_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data comes from the SeuratData package. Uncomment the InstallData()
# and install_github() lines the first time you run this; they download
# a few GB into your R library and only need to be run once.

# Setting up -------------------------------------------------------------------

# install.packages("sctransform")
# install.packages("scCustomize")
# install.packages("ggrepel")
# devtools::install_github("satijalab/seurat-data")

library(Matrix)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(scCustomize)
library(ggplot2)
library(ggrepel)

rm(list = ls())

options(future.globals.maxSize = 8000 * 1024^2)

# The data ---------------------------------------------------------------------

# SeuratData::InstallData("pbmc3k")

pbmc <- SeuratData::LoadData("pbmc3k")
pbmc

table(pbmc$seurat_annotations, useNA = "ifany")

count_mat <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")
dim(count_mat)
count_mat[11:20, 1:10]

# The two nuisances we will regress out ----------------------------------------

## Mitochondrial content
genes <- SeuratObject::Features(pbmc)
grep("^MT-", genes, value = TRUE)

pbmc <- Seurat::PercentageFeatureSet(pbmc,
                                     pattern = "^MT-",
                                     col.name = "percent.mt")
quantile(pbmc$percent.mt)

## Cell cycle
pbmc <- Seurat::NormalizeData(pbmc)

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

length(s_genes)
length(g2m_genes)
sum(s_genes %in% genes)
sum(g2m_genes %in% genes)

sum(cc.genes.updated.2019$s.genes %in% genes)
sum(cc.genes.updated.2019$g2m.genes %in% genes)

pbmc <- Seurat::CellCycleScoring(pbmc,
                                 s.features = s_genes,
                                 g2m.features = g2m_genes)
table(pbmc$Phase)

quantile(pbmc$S.Score)
quantile(pbmc$G2M.Score)

head(pbmc@meta.data)

# Running SCTransform ----------------------------------------------------------

set.seed(10)
pbmc_sct <- Seurat::SCTransform(pbmc,
                                vars.to.regress = c("percent.mt",
                                                    "S.Score",
                                                    "G2M.Score"),
                                variable.features.n = 3000,
                                seed.use = 10,
                                verbose = FALSE)
pbmc_sct

## Reading the object
SeuratObject::Assays(pbmc_sct)
SeuratObject::DefaultAssay(pbmc_sct)

SeuratObject::Layers(pbmc_sct[["SCT"]])

dim(SeuratObject::LayerData(pbmc_sct, layer = "counts", assay = "SCT"))
dim(SeuratObject::LayerData(pbmc_sct, layer = "data", assay = "SCT"))
dim(SeuratObject::LayerData(pbmc_sct, layer = "scale.data", assay = "SCT"))

data_mat <- SeuratObject::LayerData(pbmc_sct, layer = "data", assay = "SCT")
round(data_mat[1:10, 1:10], 3)

summary(Matrix::colSums(SeuratObject::LayerData(pbmc_sct, assay = "RNA",
                                                layer = "counts")))
summary(Matrix::colSums(SeuratObject::LayerData(pbmc_sct, assay = "SCT",
                                                layer = "counts")))

scale_mat <- SeuratObject::LayerData(pbmc_sct, layer = "scale.data",
                                     assay = "SCT")
round(scale_mat[1:5, 1:5], 3)
range(scale_mat)

length(Seurat::VariableFeatures(pbmc_sct))
head(Seurat::VariableFeatures(pbmc_sct), 20)

# PCA and UMAP -----------------------------------------------------------------

set.seed(10)
pbmc_sct <- Seurat::RunPCA(pbmc_sct, seed.use = 10, verbose = FALSE)
round(pbmc_sct[["pca"]]@stdev[1:20], 2)

set.seed(10)
pbmc_sct <- Seurat::RunUMAP(pbmc_sct, dims = 1:8, seed.use = 10)

## Did the regression work?
Seurat::FeaturePlot(pbmc_sct,
                    features = "percent.mt",
                    reduction = "umap")

scCustomize::FeaturePlot_scCustom(pbmc_sct,
                                  features = "percent.mt",
                                  reduction = "umap")

top_genes <- Seurat::VariableFeatures(pbmc_sct)[1:4]
top_genes

scCustomize::FeaturePlot_scCustom(pbmc_sct,
                                  features = top_genes,
                                  reduction = "umap")

Seurat::DimPlot(pbmc_sct, group.by = "seurat_annotations", label = TRUE)

# Which genes got picked, and why not the loudest ones -------------------------

count_mat <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")

count_nonzero <- Matrix::rowSums(count_mat > 0)
pct_nonzero <- (count_nonzero / ncol(count_mat)) * 100
total_count <- Matrix::rowSums(count_mat)
avg_in_nonzero <- total_count / count_nonzero

variable_genes <- Seurat::VariableFeatures(pbmc_sct)
table(variable_genes %in% rownames(count_mat))

detection_df <- data.frame(gene = rownames(count_mat),
                           avg_expression = avg_in_nonzero,
                           pct_nonzero = pct_nonzero,
                           is_variable = rownames(count_mat) %in%
                             variable_genes)

genes_to_label <- head(variable_genes, 20)

gg <- ggplot2::ggplot(detection_df,
                      ggplot2::aes(x = avg_expression, y = pct_nonzero,
                                   color = is_variable)) +
  ggplot2::geom_point(alpha = 0.5, size = 1) +
  ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                              labels = c("Other genes", "Variable features"),
                              name = "Status") +
  ggrepel::geom_text_repel(data = subset(detection_df,
                                         gene %in% genes_to_label),
                           ggplot2::aes(label = gene),
                           color = "darkred",
                           size = 3.5,
                           max.overlaps = Inf,
                           box.padding = 0.5) +
  ggplot2::labs(x = "Average expression (non-zero cells)",
                y = "Non-zero cell percentage (%)",
                title = "Expression vs. detection rate") +
  ggplot2::theme_classic()
plot(gg)

gg_zoom <- gg + ggplot2::xlim(0, 5) + ggplot2::ylim(0, 5)
plot(gg_zoom)

## Mean versus variance
n <- ncol(count_mat)
row_means <- Matrix::rowMeans(count_mat)
row_vars <- (Matrix::rowSums(count_mat^2) - n * row_means^2) / (n - 1)

head(sort(row_vars, decreasing = TRUE), 10)

variance_df <- data.frame(gene = rownames(count_mat),
                          mean = row_means,
                          variance = row_vars,
                          is_variable = rownames(count_mat) %in%
                            variable_genes)

genes_to_label <- names(sort(row_vars, decreasing = TRUE))[1:10]

gg <- ggplot2::ggplot(variance_df,
                      ggplot2::aes(x = mean, y = variance,
                                   color = is_variable)) +
  ggplot2::geom_point(alpha = 0.5, size = 1) +
  ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                              labels = c("Other genes", "Variable features"),
                              name = "Status") +
  ggrepel::geom_text_repel(data = subset(variance_df,
                                         gene %in% genes_to_label),
                           ggplot2::aes(label = gene),
                           color = "darkred",
                           size = 3.5,
                           max.overlaps = Inf,
                           box.padding = 0.5) +
  ggplot2::labs(x = "Mean expression",
                y = "Variance of expression",
                title = "Mean vs. variance") +
  ggplot2::theme_classic()
plot(gg)

gg_zoom <- gg + ggplot2::xlim(0, 5) + ggplot2::ylim(0, 20)
plot(gg_zoom)

## How different is this from the Week 1 route?
pbmc <- Seurat::FindVariableFeatures(pbmc,
                                     selection.method = "vst",
                                     nfeatures = 3000)
variable_genes_vst <- Seurat::VariableFeatures(pbmc)

length(intersect(variable_genes, variable_genes_vst))
head(variable_genes_vst, 20)

# Forcing specific genes in ----------------------------------------------------

cell_cycling <- intersect(c(s_genes, g2m_genes), SeuratObject::Features(pbmc))
length(cell_cycling)
table(cell_cycling %in% Seurat::VariableFeatures(pbmc_sct))

highly_variable_genes <- unique(c(Seurat::VariableFeatures(pbmc_sct),
                                  cell_cycling))
length(highly_variable_genes)

set.seed(10)
pbmc_final <- Seurat::SCTransform(pbmc,
                                  vars.to.regress = c("percent.mt",
                                                      "S.Score",
                                                      "G2M.Score"),
                                  residual.features = highly_variable_genes,
                                  variable.features.n = NULL,
                                  min_cells = 0,
                                  seed.use = 10,
                                  verbose = FALSE)

length(Seurat::VariableFeatures(pbmc_final))
table(cell_cycling %in% Seurat::VariableFeatures(pbmc_final))
