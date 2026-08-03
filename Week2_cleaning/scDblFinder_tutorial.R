# Week 2: Detecting doublets with scDblFinder
# Kevin Z. Lin, 2026-08-03
#
# This is the code from scDblFinder_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data downloads itself into a temporary folder, so this script should
# run as-is. Nothing is written to your machine outside that folder.

# Setting up -------------------------------------------------------------------

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("scDblFinder")

library(Seurat)
library(scDblFinder)
library(ggplot2)

rm(list = ls())

# Getting the data -------------------------------------------------------------

tmp_dir <- tempdir(check = TRUE)
base_url <- paste0("https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/",
                   "sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k")

for (component in c("filtered_feature_bc_matrix", "analysis")) {
  destination <- file.path(tmp_dir, paste0(component, ".tar.gz"))
  download.file(paste0(base_url, "_", component, ".tar.gz"),
                destfile = destination)
  untar(destination, exdir = tmp_dir)
}

## Borrowing CellRanger's clusters and UMAP
clustering_df <- read.csv(file.path(tmp_dir, "analysis", "clustering",
                                    "graphclust", "clusters.csv"),
                          row.names = 1)
umap_df <- read.csv(file.path(tmp_dir, "analysis", "umap",
                              "2_components", "projection.csv"),
                    row.names = 1)
umap_df <- as.matrix(umap_df)

head(clustering_df)
head(umap_df)

## Reading the counts
pbmc_data <- Seurat::Read10X(file.path(tmp_dir, "filtered_feature_bc_matrix"))
class(pbmc_data)
names(pbmc_data)

pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data[["Gene Expression"]],
                                   project = "pbmc1k",
                                   meta.data = clustering_df)
pbmc
head(pbmc@meta.data)

pbmc[["umap"]] <- Seurat::CreateDimReducObject(umap_df,
                                               key = "UMAP_",
                                               assay = "RNA")
pbmc

# Handing the data to Bioconductor ---------------------------------------------

Seurat::DefaultAssay(pbmc) <- "RNA"

rna_counts <- Seurat::GetAssayData(pbmc[["RNA"]], layer = "counts")
pbmc[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_counts)

sce <- Seurat::as.SingleCellExperiment(pbmc)
sce

# Running scDblFinder ----------------------------------------------------------

set.seed(10)
sce_dbl <- scDblFinder::scDblFinder(sce, clusters = NULL)

# Reading the output -----------------------------------------------------------

head(sce_dbl$scDblFinder.score)
quantile(sce_dbl$scDblFinder.score)

hist(sce_dbl$scDblFinder.score,
     breaks = 50,
     main = "Doublet scores",
     xlab = "scDblFinder score")

table(sce_dbl$scDblFinder.class)

# Back to Seurat ---------------------------------------------------------------

pbmc$scDblFinder.score <- sce_dbl$scDblFinder.score
pbmc$scDblFinder.class <- sce_dbl$scDblFinder.class

head(pbmc@meta.data)

all(colnames(pbmc) == colnames(sce_dbl))

# Looking at the doublets ------------------------------------------------------

Seurat::DimPlot(pbmc, group.by = "scDblFinder.class") +
  ggplot2::ggtitle("Doublet calls")

Seurat::FeaturePlot(pbmc, features = "scDblFinder.score") +
  ggplot2::ggtitle("Doublet scores")

Seurat::DimPlot(pbmc, group.by = "Cluster", label = TRUE) +
  ggplot2::ggtitle("CellRanger clusters")

# Filtering, and the result ----------------------------------------------------

pbmc <- subset(pbmc, scDblFinder.class == "singlet")
pbmc

Seurat::DimPlot(pbmc, group.by = "Cluster", label = TRUE) +
  ggplot2::ggtitle("Singlets only")
