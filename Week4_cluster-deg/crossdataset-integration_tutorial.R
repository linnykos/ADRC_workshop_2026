# Week 4: Batch correction, part 2 — integrating two studies, and correcting the expression matrix
# Kevin Z. Lin, 2026-08-17
#
# This is the code from crossdataset-integration_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# Two datasets (seaad_microglia.RData, 525 MB, and rosmap_microglia.RData,
# 1.5 GB) download themselves into a temporary folder. Nothing is written to
# your machine outside that folder. The first run takes several minutes on
# the downloads.
#
# UNLIKE THE OTHER SCRIPTS HERE, THIS ONE DOES NOT RUN AS-IS. The two
# integration calls and the label transfer are meant to run on a server, via
# crossdataset-integration_server.R and its .slurm file; here they are left
# in place so you can read them, but running them on a laptop takes upwards
# of half an hour and a lot of memory. The load() of
# crossdataset-integration_results.RData part-way down is where the server's
# output comes back in. See the .html for the whole story.

# Setting up -------------------------------------------------------------------

# install.packages(c("Seurat", "RANN", "ggplot2", "patchwork"))

library(Seurat)
library(SeuratObject)
library(RANN)
library(ggplot2)
library(patchwork)

rm(list = ls())

# The data: two studies, one brain region --------------------------------------

data_dir <- tempdir(check = TRUE)

# If you have already downloaded either file, point data_dir at the folder
# holding them and the downloads below will be skipped, e.g.:
# data_dir <- "/path/to/your/data"

options(timeout = 3600)

seaad_url <- paste0("https://www.dropbox.com/scl/fi/cei02rwh3gx61uu9zzeql/",
                    "seaad_microglia.RData",
                    "?rlkey=04erfwn3w2lrnr80pefwfh8fu&dl=1")
seaad_file <- file.path(data_dir, "seaad_microglia.RData")
seaad_size <- 525672784

if (!file.exists(seaad_file) || file.size(seaad_file) != seaad_size) {
  download.file(seaad_url, destfile = seaad_file, mode = "wb")
}

file.size(seaad_file) == seaad_size

rosmap_url <- paste0("https://www.dropbox.com/scl/fi/siecgiljoks48pkin42ts/",
                     "rosmap_microglia.RData",
                     "?rlkey=vwnmoeynhbsuxkiumz36hy301&dl=1")
rosmap_file <- file.path(data_dir, "rosmap_microglia.RData")
rosmap_size <- 1507998516

if (!file.exists(rosmap_file) || file.size(rosmap_file) != rosmap_size) {
  download.file(rosmap_url, destfile = rosmap_file, mode = "wb")
}

file.size(rosmap_file) == rosmap_size

load(seaad_file)
load(rosmap_file)

ls()

seaad <- seurat_obj
rosmap <- rosmap_obj
rm(list = c("seurat_obj", "rosmap_obj")); gc(TRUE)

seaad
rosmap

SeuratObject::Layers(seaad[["RNA"]])
length(unique(seaad$donor_id))
table(seaad$Supertype)

SeuratObject::Layers(rosmap[["RNA"]])
length(unique(rosmap$donor_id))
table(rosmap$state)
table(rosmap$sequencing_batch)
table(rosmap$ad_diagnosis)

# Making the two objects comparable --------------------------------------------

## Barcodes
utils::head(Seurat::Cells(seaad), 3)
utils::head(Seurat::Cells(rosmap), 3)

seaad <- Seurat::RenameCells(seaad, add.cell.id = "seaad")
rosmap <- Seurat::RenameCells(rosmap, add.cell.id = "rosmap")

utils::head(Seurat::Cells(seaad), 3)
utils::head(Seurat::Cells(rosmap), 3)

seaad$dataset <- "SEA-AD"
rosmap$dataset <- "ROSMAP"

## Genes
length(rownames(seaad))
length(rownames(rosmap))

common_genes <- intersect(rownames(seaad), rownames(rosmap))
length(common_genes)

length(setdiff(rownames(seaad), common_genes))
length(setdiff(rownames(rosmap), common_genes))

seaad <- subset(seaad, features = common_genes)
rosmap <- subset(rosmap, features = common_genes)

dim(seaad)
dim(rosmap)

## Size
ncol(seaad)
ncol(rosmap)

set.seed(10)
seaad_cells <- sample(Seurat::Cells(seaad), size = 10000)
seaad <- subset(seaad, cells = seaad_cells)

ncol(seaad)
ncol(rosmap)

## Normalization and variable genes, within each study
seaad <- Seurat::NormalizeData(seaad, verbose = FALSE)
rosmap <- Seurat::NormalizeData(rosmap, verbose = FALSE)

seaad <- Seurat::FindVariableFeatures(seaad, selection.method = "vst",
                                      nfeatures = 2000, verbose = FALSE)
rosmap <- Seurat::FindVariableFeatures(rosmap, selection.method = "vst",
                                       nfeatures = 2000, verbose = FALSE)

length(intersect(Seurat::VariableFeatures(seaad),
                 Seurat::VariableFeatures(rosmap)))

integration_features <- Seurat::SelectIntegrationFeatures(
  object.list = list(seaad, rosmap),
  nfeatures = 2000,
  verbose = FALSE)

length(integration_features)
utils::head(integration_features, 20)

# What the two studies look like uncorrected -----------------------------------

merged <- merge(seaad, rosmap)
merged

SeuratObject::Layers(merged[["RNA"]])
table(merged$dataset)

merged <- Seurat::NormalizeData(merged, verbose = FALSE)
merged <- Seurat::FindVariableFeatures(merged, selection.method = "vst",
                                       nfeatures = 2000, verbose = FALSE)
merged <- Seurat::ScaleData(merged, verbose = FALSE)

set.seed(10)
merged <- Seurat::RunPCA(merged,
                         features = Seurat::VariableFeatures(merged),
                         npcs = 30,
                         reduction.name = "pca.uncorrected",
                         verbose = FALSE)

set.seed(10)
merged <- Seurat::RunUMAP(merged,
                          reduction = "pca.uncorrected",
                          dims = 1:30,
                          reduction.name = "umap.uncorrected",
                          seed.use = 10,
                          verbose = FALSE)

Seurat::Reductions(merged)

gg <- Seurat::DimPlot(merged, reduction = "umap.uncorrected",
                      group.by = "dataset", raster = TRUE)
gg <- gg + ggplot2::ggtitle("Uncorrected, coloured by study")
plot(gg)

gg <- Seurat::DimPlot(merged, reduction = "umap.uncorrected",
                      group.by = "dataset", split.by = "dataset",
                      raster = TRUE)
plot(gg)

## Measuring the separation instead of describing it
neighbor_mixing <- function(embedding_mat, batch_vec, k = 30){
  # nn2 returns each point as its own first neighbour, so drop column 1.
  nn_out <- RANN::nn2(data = embedding_mat, k = k + 1)
  nn_idx <- nn_out$nn.idx[, -1, drop = FALSE]

  # batch_vec[nn_idx] unwinds column-major, so rebuild the matrix with the same
  # nrow to keep cell i on row i.
  neighbor_batch_mat <- matrix(batch_vec[nn_idx], nrow = nrow(nn_idx))

  rowMeans(neighbor_batch_mat != batch_vec)
}

dataset_vec <- merged$dataset
dataset_share <- table(dataset_vec) / length(dataset_vec)
dataset_share

# Under perfect mixing, a cell's neighbours are a random sample of everything
# else, so the expected other-study fraction is one minus its own study's share.
1 - dataset_share

mixing_uncorrected <- neighbor_mixing(
  embedding_mat = Seurat::Embeddings(merged, "pca.uncorrected")[, 1:30],
  batch_vec = dataset_vec,
  k = 30)

tapply(mixing_uncorrected, dataset_vec, stats::median)

# Integrating so that the expression matrix is corrected too -------------------

# Runs on the server, not here. See crossdataset-integration_server.R.
integration_anchors <- Seurat::FindIntegrationAnchors(
  object.list = list(seaad, rosmap),
  anchor.features = integration_features,
  reduction = "cca",
  dims = 1:30,
  verbose = FALSE)

class(integration_anchors)

# Runs on the server, not here. See crossdataset-integration_server.R.
integrated <- Seurat::IntegrateData(anchorset = integration_anchors,
                                    new.assay.name = "integrated",
                                    dims = 1:30,
                                    verbose = FALSE)

# Bringing the result back -----------------------------------------------------

results_file <- "crossdataset-integration_results.RData"

# The default looks for the file next to this document. If you put it somewhere
# else, point at it here, e.g.:
# results_file <- "/path/to/crossdataset-integration_results.RData"

if (!file.exists(results_file)) {
  stop("Cannot find ", results_file, ". Run crossdataset-integration_server.R ",
       "on a server (see the previous section) and copy the file it writes ",
       "next to this document.", call. = FALSE)
}

local_seaad_cells <- seaad_cells
local_features <- integration_features

load(results_file)

ls()
run_info

# Did both machines draw the same 10,000 nuclei, and select the same genes?
identical(local_seaad_cells, seaad_cells)
identical(local_features, integration_features)

# And does the object that came back describe the cells we have in front of us?
setequal(Seurat::Cells(integrated),
         c(Seurat::Cells(seaad), Seurat::Cells(rosmap)))

stopifnot(identical(local_seaad_cells, seaad_cells),
          setequal(Seurat::Cells(integrated),
                   c(Seurat::Cells(seaad), Seurat::Cells(rosmap))))

integrated
Seurat::Assays(integrated)
Seurat::DefaultAssay(integrated)

SeuratObject::Layers(integrated[["RNA"]])

integrated[["RNA"]] <- SeuratObject::JoinLayers(integrated[["RNA"]])

SeuratObject::Layers(integrated[["RNA"]])

# What is actually inside the `integrated` assay -------------------------------

class(integrated[["integrated"]])
SeuratObject::Layers(integrated[["integrated"]])

dim(integrated[["RNA"]])
dim(integrated[["integrated"]])

integrated_data <- SeuratObject::LayerData(integrated[["integrated"]],
                                           layer = "data")
rna_data <- SeuratObject::LayerData(integrated[["RNA"]], layer = "data")

stopifnot(all(colnames(integrated_data) == colnames(rna_data)))

range(rna_data)
range(integrated_data)
sum(integrated_data < 0) / length(integrated_data)

marker_candidates <- c("P2RY12", "CX3CR1", "CSF1R", "C1QB", "APOE", "SPP1")
marker_gene <- marker_candidates[marker_candidates %in%
                                   rownames(integrated_data)][1]
marker_gene

sum(rna_data[marker_gene, ] == 0)
sum(integrated_data[marker_gene, ] == 0)

compare_df <- data.frame(
  rna = rna_data[marker_gene, ],
  integrated = integrated_data[marker_gene, ],
  dataset = integrated$dataset)

gg <- ggplot2::ggplot(compare_df,
                      ggplot2::aes(x = rna, y = integrated, color = dataset))
gg <- gg + ggplot2::geom_point(size = 0.3, alpha = 0.2)
gg <- gg + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed")
gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "red")
gg <- gg + ggplot2::labs(x = paste0(marker_gene, ", measured (log-normalized)"),
                         y = paste0(marker_gene, ", after integration"),
                         title = "The same gene, the same cells, two matrices")
plot(gg)

# Re-embedding on the corrected values -----------------------------------------

Seurat::DefaultAssay(integrated) <- "integrated"

integrated <- Seurat::ScaleData(integrated, verbose = FALSE)

set.seed(10)
integrated <- Seurat::RunPCA(integrated,
                             npcs = 30,
                             reduction.name = "pca.integrated",
                             verbose = FALSE)

set.seed(10)
integrated <- Seurat::RunUMAP(integrated,
                              reduction = "pca.integrated",
                              dims = 1:30,
                              reduction.name = "umap.integrated",
                              seed.use = 10,
                              verbose = FALSE)

Seurat::Reductions(integrated)

gg1 <- Seurat::DimPlot(merged, reduction = "umap.uncorrected",
                       group.by = "dataset", raster = TRUE) +
  ggplot2::ggtitle("Before")
gg2 <- Seurat::DimPlot(integrated, reduction = "umap.integrated",
                       group.by = "dataset", raster = TRUE) +
  ggplot2::ggtitle("After")
plot(gg1 + gg2)

mixing_integrated <- neighbor_mixing(
  embedding_mat = Seurat::Embeddings(integrated, "pca.integrated")[, 1:30],
  batch_vec = integrated$dataset,
  k = 30)

tapply(mixing_uncorrected, merged$dataset, stats::median)
tapply(mixing_integrated, integrated$dataset, stats::median)
1 - dataset_share

mixing_df <- rbind(
  data.frame(mixing = mixing_uncorrected, dataset = merged$dataset,
             stage = "before"),
  data.frame(mixing = mixing_integrated, dataset = integrated$dataset,
             stage = "after"))
mixing_df$stage <- factor(mixing_df$stage, levels = c("before", "after"))

gg <- ggplot2::ggplot(mixing_df,
                      ggplot2::aes(x = mixing, fill = dataset))
gg <- gg + ggplot2::geom_histogram(bins = 40, alpha = 0.6,
                                   position = "identity")
gg <- gg + ggplot2::facet_wrap(~ stage)
gg <- gg + ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed")
gg <- gg + ggplot2::labs(x = "fraction of 30 neighbours from the other study",
                         y = "cells")
plot(gg)

# Going back to the counts -----------------------------------------------------

Seurat::DefaultAssay(integrated) <- "RNA"
Seurat::DefaultAssay(integrated)

SeuratObject::Layers(integrated[["RNA"]])

Seurat::DefaultAssay(integrated) <- "RNA"
gg1 <- Seurat::FeaturePlot(integrated, features = marker_gene,
                           reduction = "umap.integrated", raster = TRUE) +
  ggplot2::ggtitle(paste0(marker_gene, ", measured"))

Seurat::DefaultAssay(integrated) <- "integrated"
gg2 <- Seurat::FeaturePlot(integrated, features = marker_gene,
                           reduction = "umap.integrated", raster = TRUE) +
  ggplot2::ggtitle(paste0(marker_gene, ", corrected"))

Seurat::DefaultAssay(integrated) <- "RNA"
plot(gg1 + gg2)

# The alternative that never edits expression: label transfer ------------------

# Runs on the server, not here. See crossdataset-integration_server.R.
transfer_anchors <- Seurat::FindTransferAnchors(reference = seaad,
                                                query = rosmap,
                                                dims = 1:30,
                                                verbose = FALSE)

predictions_df <- Seurat::TransferData(anchorset = transfer_anchors,
                                       refdata = seaad$Supertype,
                                       dims = 1:30,
                                       verbose = FALSE)

dim(predictions_df)
utils::head(colnames(predictions_df))
utils::head(predictions_df[, c("predicted.id", "prediction.score.max")], 5)

rosmap <- Seurat::AddMetaData(rosmap, metadata = predictions_df)

table(rosmap$predicted.id)
stats::quantile(rosmap$prediction.score.max)

label_tab <- table(rosmap$state, rosmap$predicted.id)
label_tab

round(label_tab / rowSums(label_tab), 2)
