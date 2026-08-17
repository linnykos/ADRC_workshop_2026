# Week 4: Clustering, cluster stability, and marker genes
# Kevin Z. Lin, 2026-08-15
#
# This is the code from clustering-markers_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data (seaad_microglia.RData, 525 MB) downloads itself into a temporary
# folder, so this script should run as-is. Nothing is written to your machine
# outside that folder. The first run takes a few minutes on the download.

# Setting up -------------------------------------------------------------------

# install.packages(c("Seurat", "ggplot2", "pheatmap", "patchwork"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("bluster", "EnhancedVolcano"))
# devtools::install_github("immunogenomics/presto")  # makes wilcox much faster

library(Seurat)
library(SeuratObject)
library(bluster)
library(EnhancedVolcano)
library(pheatmap)
library(patchwork)
library(ggplot2)

rm(list = ls())

# The data ---------------------------------------------------------------------

data_dir <- tempdir(check = TRUE)

# If you have already downloaded seaad_microglia.RData, point data_dir at the
# folder holding it instead and the download below will be skipped, e.g.:
# data_dir <- "/path/to/your/data"

options(timeout = 3600)

data_url <- paste0("https://www.dropbox.com/scl/fi/cei02rwh3gx61uu9zzeql/",
                   "seaad_microglia.RData",
                   "?rlkey=04erfwn3w2lrnr80pefwfh8fu&dl=1")
data_file <- file.path(data_dir, "seaad_microglia.RData")
expected_size <- 525672784

if (!file.exists(data_file) || file.size(data_file) != expected_size) {
  download.file(data_url, destfile = data_file, mode = "wb")
}

file.size(data_file) == expected_size

load(data_file)
seurat_obj

SeuratObject::Layers(seurat_obj[["RNA"]])
Seurat::Reductions(seurat_obj)

dim(Seurat::Embeddings(seurat_obj, "scVI"))
dim(Seurat::Embeddings(seurat_obj, "umap"))

length(unique(seurat_obj$donor_id))
table(seurat_obj$assay)
table(seurat_obj$Supertype)

gg1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Supertype",
                       raster = TRUE) +
  ggplot2::ggtitle("SEA-AD supertypes")
gg2 <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "donor_id",
                       raster = TRUE) +
  Seurat::NoLegend() +
  ggplot2::ggtitle("Donor (89 of them)")
plot(gg1 + gg2)

# The naive route: cluster without correcting for donor ------------------------

seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           selection.method = "vst",
                                           nfeatures = 2000,
                                           verbose = FALSE)
seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)

set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             features = Seurat::VariableFeatures(seurat_obj),
                             npcs = 30,
                             verbose = FALSE)

seurat_obj[["RNA"]]$scale.data <- NULL

head(Seurat::VariableFeatures(seurat_obj), 20)

gg <- Seurat::ElbowPlot(seurat_obj, ndims = 30) +
  ggplot2::ggtitle("Uncorrected PCA")
plot(gg)

pca_stdev <- Seurat::Stdev(seurat_obj, reduction = "pca")
stdev_drop <- abs(diff(pca_stdev))

num_pcs <- min(which(stdev_drop < 0.1))
num_pcs

round(pca_stdev[1:15], 2)

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj,
                                    reduction = "pca",
                                    dims = 1:num_pcs,
                                    verbose = FALSE)
seurat_obj <- Seurat::FindClusters(seurat_obj,
                                   resolution = 0.5,
                                   cluster.name = "pca_clusters",
                                   verbose = FALSE)

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              reduction = "pca",
                              dims = 1:num_pcs,
                              reduction.name = "umap_pca",
                              verbose = FALSE)

table(seurat_obj$pca_clusters)

gg1 <- Seurat::DimPlot(seurat_obj, reduction = "umap_pca",
                       group.by = "pca_clusters", label = TRUE,
                       raster = TRUE) +
  Seurat::NoLegend() +
  ggplot2::ggtitle("Uncorrected: clusters")
gg2 <- Seurat::DimPlot(seurat_obj, reduction = "umap_pca",
                       group.by = "donor_id", raster = TRUE) +
  Seurat::NoLegend() +
  ggplot2::ggtitle("Uncorrected: donor")
plot(gg1 + gg2)

## Measuring the damage instead of eyeballing it
donor_purity <- function(clusters, donors) {
  cluster_levels <- sort(unique(as.character(clusters)))
  out_df <- data.frame(cluster = cluster_levels,
                       n_cells = NA_integer_,
                       top_donor_frac = NA_real_,
                       donors_for_half = NA_integer_)

  for (i in seq_along(cluster_levels)) {
    in_cluster <- as.character(clusters) == cluster_levels[i]
    donor_tab <- sort(table(donors[in_cluster]), decreasing = TRUE)
    out_df[i, "n_cells"] <- sum(in_cluster)
    out_df[i, "top_donor_frac"] <- max(donor_tab) / sum(donor_tab)
    out_df[i, "donors_for_half"] <-
      min(which(cumsum(donor_tab) >= 0.5 * sum(donor_tab)))
  }
  out_df
}

purity_pca <- donor_purity(seurat_obj$pca_clusters, seurat_obj$donor_id)
purity_pca

# The corrected route: cluster on the integrated embedding ---------------------

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj,
                                    reduction = "scVI",
                                    dims = 1:10,
                                    verbose = FALSE)
seurat_obj <- Seurat::FindClusters(seurat_obj,
                                   resolution = 0.5,
                                   cluster.name = "scvi_clusters",
                                   verbose = FALSE)

table(seurat_obj$scvi_clusters)

gg1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                       group.by = "scvi_clusters", label = TRUE,
                       raster = TRUE) +
  Seurat::NoLegend() +
  ggplot2::ggtitle("Integrated: clusters")
gg2 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                       group.by = "donor_id", raster = TRUE) +
  Seurat::NoLegend() +
  ggplot2::ggtitle("Integrated: donor")
plot(gg1 + gg2)

purity_scvi <- donor_purity(seurat_obj$scvi_clusters, seurat_obj$donor_id)
purity_scvi

summary(purity_pca[, "top_donor_frac"])
summary(purity_scvi[, "top_donor_frac"])

sum(purity_pca[, "donors_for_half"] <= 3)
sum(purity_scvi[, "donors_for_half"] <= 3)

# Do the two clusterings agree? ------------------------------------------------

concordance_tab <- table(paste0("pca_", seurat_obj$pca_clusters),
                         paste0("scvi_", seurat_obj$scvi_clusters))
concordance_tab

pheatmap::pheatmap(log10(concordance_tab + 10),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = "Uncorrected vs integrated clusters (log10 counts)")

# Are the clusters stable? -----------------------------------------------------

set.seed(10)
stability_cells <- sample(colnames(seurat_obj), 10000)
scvi_embedding <- Seurat::Embeddings(seurat_obj, "scVI")[stability_cells, ]
dim(scvi_embedding)

cluster_fun <- function(x) {
  graph_obj <- bluster::makeSNNGraph(x, type = "jaccard")
  igraph::cluster_louvain(graph_obj)$membership
}

set.seed(10)
original_clusters <- cluster_fun(scvi_embedding)
table(original_clusters)

set.seed(10)
stability_mat <- bluster::bootstrapStability(scvi_embedding,
                                             FUN = cluster_fun,
                                             clusters = original_clusters)
dim(stability_mat)
round(diag(stability_mat), 2)

pheatmap::pheatmap(stability_mat,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = "Bootstrap cluster coherence")

# Choosing a resolution --------------------------------------------------------

set.seed(10)
seurat_obj <- Seurat::FindClusters(seurat_obj,
                                   resolution = c(0.1, 0.2, 0.5, 1.0),
                                   verbose = FALSE)

resolution_cols <- c("RNA_snn_res.0.1", "RNA_snn_res.0.2",
                     "RNA_snn_res.0.5", "RNA_snn_res.1")
sapply(resolution_cols,
       function(x) length(unique(seurat_obj@meta.data[, x])))

gg_list <- lapply(resolution_cols, function(x) {
  Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = x,
                  label = TRUE, raster = TRUE) +
    Seurat::NoLegend() +
    ggplot2::ggtitle(x)
})
plot(patchwork::wrap_plots(gg_list, ncol = 2))

## Where are the cells actually dense?
umap_df <- Seurat::FetchData(seurat_obj, vars = c("umap_1", "umap_2"))

gg1 <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2)) +
  ggplot2::geom_point(alpha = 0.15, size = 0.01) +
  ggplot2::geom_density2d(colour = "red", bins = 12) +
  ggplot2::ggtitle("Cell density, contour lines") +
  ggplot2::theme_bw()

gg2 <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2)) +
  ggplot2::stat_density_2d(geom = "polygon", contour = TRUE,
                           ggplot2::aes(fill = ggplot2::after_stat(level)),
                           bins = 12, alpha = 0.4) +
  ggplot2::scale_fill_distiller(palette = "Blues", direction = 1) +
  ggplot2::geom_density2d(colour = "black", bins = 12, linewidth = 0.2) +
  ggplot2::ggtitle("Cell density, filled contours") +
  ggplot2::theme_bw()

plot(gg1 + gg2)

set.seed(10)
donor_subset <- sample(unique(seurat_obj$donor_id), 6)

donor_df <- Seurat::FetchData(seurat_obj,
                              vars = c("umap_1", "umap_2", "donor_id"))
donor_df <- donor_df[donor_df$donor_id %in% donor_subset, ]

gg <- ggplot2::ggplot(donor_df, ggplot2::aes(x = umap_1, y = umap_2)) +
  ggplot2::geom_point(data = umap_df, colour = "grey85", size = 0.01) +
  ggplot2::geom_density2d(colour = "red", bins = 8) +
  ggplot2::facet_wrap(~ donor_id, ncol = 3) +
  ggplot2::ggtitle("Density per donor, six donors, whole dataset in grey") +
  ggplot2::theme_bw()
plot(gg)

purity_by_resolution <- lapply(resolution_cols, function(x) {
  purity_df <- donor_purity(seurat_obj@meta.data[, x], seurat_obj$donor_id)
  data.frame(resolution = x,
             n_clusters = nrow(purity_df),
             median_top_donor = median(purity_df[, "top_donor_frac"]),
             max_top_donor = max(purity_df[, "top_donor_frac"]),
             min_cluster_size = min(purity_df[, "n_cells"]))
})
do.call(rbind, purity_by_resolution)

cells_tab <- table(seurat_obj$RNA_snn_res.0.2, seurat_obj$donor_id)
apply(cells_tab, 1,
      function(x) c(donors_with_10_or_more = sum(x >= 10),
                    median_cells = median(x)))

Seurat::Idents(seurat_obj) <- "RNA_snn_res.0.2"
table(Seurat::Idents(seurat_obj))

# Comparing two clusters -------------------------------------------------------

de_res <- Seurat::FindMarkers(seurat_obj,
                              ident.1 = "0",
                              ident.2 = "1",
                              test.use = "wilcox")
head(de_res, 10)
nrow(de_res)

bh_adj <- p.adjust(de_res$p_val, method = "BH")
sum(bh_adj < 0.05)

p_threshold <- max(de_res$p_val[bh_adj < 0.05])
p_threshold

gg <- EnhancedVolcano::EnhancedVolcano(
  de_res,
  lab = rownames(de_res),
  x = "avg_log2FC",
  y = "p_val",
  pCutoff = p_threshold,
  FCcutoff = 0.5,
  title = "Cluster 0 versus cluster 1",
  subtitle = "Wilcoxon, cells as units; line at BH-adjusted p < 0.05")
plot(gg)

# Marker genes -----------------------------------------------------------------

all_markers <- Seurat::FindAllMarkers(seurat_obj,
                                      only.pos = TRUE,
                                      test.use = "wilcox",
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25,
                                      verbose = FALSE)
dim(all_markers)
table(all_markers$cluster)

top_markers <- do.call(rbind, lapply(split(all_markers, all_markers$cluster),
                                     function(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 5)
}))
top_markers[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

marker_genes <- unique(top_markers[, "gene"])
gg <- Seurat::DotPlot(seurat_obj, features = marker_genes) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                     vjust = 0.5)) +
  ggplot2::ggtitle("Top 5 markers per cluster by fold change")
plot(gg)

## Ranking by fold change is not the same as ranking by p-value
cluster0_markers <- all_markers[all_markers$cluster == "0", ]

by_pval <- head(cluster0_markers[order(cluster0_markers$p_val), "gene"], 10)
by_fc <- head(cluster0_markers[order(cluster0_markers$avg_log2FC,
                                     decreasing = TRUE), "gene"], 10)

by_pval
by_fc
length(intersect(by_pval, by_fc))

## A curated panel beats a derived list for naming
panel_list <- list(
  homeostatic = c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", "CSF1R", "SALL1"),
  activated = c("APOE", "SPP1", "GPNMB", "LPL", "ITGAX", "CD9"),
  antigen = c("CD74", "HLA-DRA", "HLA-DRB1"),
  pvm = c("CD163", "MRC1", "F13A1", "LYVE1"),
  lymphocyte = c("CD2", "CD3E", "IL7R", "SKAP1", "THEMIS", "CCL5"),
  proliferating = c("MKI67", "TOP2A"),
  other_lineage = c("SNAP25", "PLP1", "SLC1A2", "CLDN5")
)

panel_genes <- unlist(panel_list)
setdiff(panel_genes, rownames(seurat_obj))
panel_genes <- intersect(panel_genes, rownames(seurat_obj))

gg <- Seurat::DotPlot(seurat_obj, features = panel_genes) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                     vjust = 0.5)) +
  ggplot2::ggtitle("Curated microglial panels")
plot(gg)

gg <- Seurat::FeaturePlot(seurat_obj, reduction = "umap",
                          features = c("P2RY12", "CX3CR1", "APOE", "SPP1",
                                       "CD163", "MRC1"),
                          raster = TRUE, ncol = 3)
plot(gg)

# Checking against a reference annotation --------------------------------------

supertype_tab <- table(paste0("cluster_", Seurat::Idents(seurat_obj)),
                       seurat_obj$Supertype)
supertype_tab

pheatmap::pheatmap(log10(supertype_tab + 10),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = "Our clusters vs SEA-AD supertypes (log10 counts)")

round(prop.table(supertype_tab, margin = 1), 2)
