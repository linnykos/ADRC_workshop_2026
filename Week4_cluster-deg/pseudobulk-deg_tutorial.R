# Week 4: Differential expression across donors, and what the genes mean
# Kevin Z. Lin, 2026-08-15
#
# This is the code from pseudobulk-deg_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data (seaad_microglia.RData, 525 MB) downloads itself into a temporary
# folder, so this script should run as-is. Nothing is written to your machine
# outside that folder. The first run takes a few minutes on the download.

# Setting up -------------------------------------------------------------------

# install.packages(c("Seurat", "ggplot2", "pheatmap", "RColorBrewer"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano",
#                        "clusterProfiler", "org.Hs.eg.db", "enrichplot",
#                        "variancePartition"))

library(Seurat)
library(SeuratObject)
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(variancePartition)
library(ggplot2)
library(pheatmap)

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

ls()

head(seurat_obj@meta.data[, c("donor_id", "Supertype", "assay", "sex",
                              "ADNC", "APOE4.status", "Cognitive.status")])

table(seurat_obj$ADNC)
table(seurat_obj$assay)

## Two subsets, both deliberate
seurat_obj <- subset(seurat_obj,
                     ADNC %in% c("High", "Intermediate", "Low", "Not AD"))
seurat_obj <- subset(seurat_obj, assay == "10x 3' v3")
seurat_obj

length(unique(seurat_obj$donor_id))

# Choosing a population to test ------------------------------------------------

cells_tab <- table(seurat_obj$Supertype, seurat_obj$donor_id)
dim(cells_tab)

apply(cells_tab, 1,
      function(x) c(donors_with_10_or_more = sum(x >= 10),
                    median_cells = median(x),
                    max_cells = max(x)))

seurat_obj <- subset(seurat_obj, Supertype == "Micro-PVM_2")
seurat_obj

summary(as.vector(table(seurat_obj$donor_id)))

# Aggregating to pseudobulk ----------------------------------------------------

pseudo_seurat <- Seurat::AggregateExpression(
  seurat_obj,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("donor_id", "APOE4.status", "ADNC", "sex"),
  verbose = FALSE)

pseudo_seurat
head(pseudo_seurat@meta.data)

mat_pseudobulk <- SeuratObject::LayerData(pseudo_seurat,
                                          layer = "counts",
                                          assay = "RNA")
metadata_pseudobulk <- pseudo_seurat@meta.data

dim(mat_pseudobulk)
mat_pseudobulk[1:5, 1:3]

## How many cells is each column made of?
cells_per_donor <- table(seurat_obj$donor_id)
metadata_pseudobulk[, "n_cells"] <-
  as.vector(cells_per_donor[metadata_pseudobulk[, "donor_id"]])

summary(metadata_pseudobulk[, "n_cells"])
sort(metadata_pseudobulk[, "n_cells"])[1:8]

min_cells <- 10
keep_samples <- metadata_pseudobulk[, "n_cells"] >= min_cells
sum(!keep_samples)

mat_pseudobulk <- mat_pseudobulk[, keep_samples]
metadata_pseudobulk <- metadata_pseudobulk[keep_samples, ]
dim(mat_pseudobulk)

## Preparing the variables
metadata_pseudobulk[, "sex"] <- factor(metadata_pseudobulk[, "sex"])
metadata_pseudobulk[, "APOE4.status"] <-
  factor(metadata_pseudobulk[, "APOE4.status"])

adnc_vec <- rep("NonAD", nrow(metadata_pseudobulk))
adnc_vec[metadata_pseudobulk[, "ADNC"] %in% c("High", "Intermediate")] <- "AD"
metadata_pseudobulk[, "ADNC"] <- relevel(factor(adnc_vec), ref = "NonAD")

summary(metadata_pseudobulk[, c("sex", "APOE4.status", "ADNC")])

## Check the alignment before you trust it
all(colnames(mat_pseudobulk) == rownames(metadata_pseudobulk))

# Choosing which genes to test -------------------------------------------------

min_samples <- ceiling(0.25 * ncol(mat_pseudobulk))
min_samples

genes_detected <- Matrix::rowSums(mat_pseudobulk >= 5) >= min_samples
sum(genes_detected)

mat_pseudobulk <- mat_pseudobulk[genes_detected, ]
dim(mat_pseudobulk)

# The design -------------------------------------------------------------------

design_formula <- ~ sex + APOE4.status + ADNC
design_mat <- model.matrix(design_formula, metadata_pseudobulk)

colnames(design_mat)
dim(design_mat)

qr(design_mat)$rank
ncol(design_mat)
qr(design_mat)$rank == ncol(design_mat)

table(metadata_pseudobulk[, "ADNC"], metadata_pseudobulk[, "sex"])
table(metadata_pseudobulk[, "ADNC"], metadata_pseudobulk[, "APOE4.status"])

# Looking at the data before testing it ----------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = design_formula)
dds

## Why the counts need transforming first
vsd <- DESeq2::vst(dds, blind = TRUE)

## What to look for in the correlation heatmap
vsd_mat <- SummarizedExperiment::assay(vsd)
cor_mat <- cor(vsd_mat)

annotation_df <- data.frame(ADNC = metadata_pseudobulk[, "ADNC"],
                            sex = metadata_pseudobulk[, "sex"],
                            row.names = rownames(metadata_pseudobulk))

pheatmap::pheatmap(cor_mat,
                   annotation_col = annotation_df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Pseudobulk sample-sample correlation")

cor_offdiag <- cor_mat[lower.tri(cor_mat)]
round(quantile(cor_offdiag, c(0, 0.05, 0.5, 0.95, 1)), 3)

mean_cor <- (rowSums(cor_mat) - 1) / (ncol(cor_mat) - 1)
round(sort(mean_cor)[1:5], 3)

metadata_pseudobulk[names(sort(mean_cor)[1:5]), "n_cells"]

cor(mean_cor, log10(metadata_pseudobulk[, "n_cells"]))

## The same questions, in a PCA
gg <- DESeq2::plotPCA(vsd, intgroup = "ADNC") +
  ggplot2::ggtitle("Pseudobulk PCA, coloured by ADNC")
plot(gg)

pca_df <- DESeq2::plotPCA(vsd, intgroup = "ADNC", returnData = TRUE)
pca_df[, "n_cells"] <- metadata_pseudobulk[, "n_cells"]

gg <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_point(ggplot2::aes(colour = n_cells), size = 3) +
  ggplot2::scale_colour_viridis_c(trans = "log10") +
  ggplot2::ggtitle("The same PCA, coloured by cells per donor")
plot(gg)

cor(pca_df[, "PC1"], log10(pca_df[, "n_cells"]))
cor(pca_df[, "PC2"], log10(pca_df[, "n_cells"]))

summary(aov(pca_df[, "PC1"] ~ metadata_pseudobulk[, "ADNC"]))

# Testing ----------------------------------------------------------------------

dds <- DESeq2::DESeq(dds)
DESeq2::resultsNames(dds)

deseq2_res <- DESeq2::results(dds, name = "ADNC_AD_vs_NonAD")
summary(deseq2_res)

deseq2_res_ordered <- deseq2_res[order(deseq2_res$padj), ]
head(as.data.frame(deseq2_res_ordered), 10)

## Shrinking the fold changes
deseq2_shrunk <- DESeq2::lfcShrink(dds,
                                   coef = "ADNC_AD_vs_NonAD",
                                   type = "apeglm")
head(as.data.frame(deseq2_shrunk[order(deseq2_shrunk$padj), ]), 10)

par(mfrow = c(1, 2))
DESeq2::plotMA(deseq2_res, ylim = c(-3, 3), main = "Unshrunk")
DESeq2::plotMA(deseq2_shrunk, ylim = c(-3, 3), main = "Shrunk (apeglm)")
par(mfrow = c(1, 1))

gg <- EnhancedVolcano::EnhancedVolcano(
  deseq2_shrunk,
  lab = rownames(deseq2_shrunk),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  title = "Micro-PVM_2, AD vs non-AD",
  subtitle = "Pseudobulk by donor, adjusted for sex and APOE4 status")
plot(gg)

# Interrogating the answer -----------------------------------------------------

top_genes <- rownames(deseq2_res_ordered)[1:20]
top_genes

sum(grepl("^MT-", top_genes))

mito_genes <- grep("^MT-", rownames(mat_pseudobulk), value = TRUE)
mito_genes

mito_fraction <- Matrix::colSums(mat_pseudobulk[mito_genes, ]) /
  Matrix::colSums(mat_pseudobulk)
metadata_pseudobulk[, "mito_fraction"] <- mito_fraction

summary(mito_fraction)
tapply(mito_fraction, metadata_pseudobulk[, "ADNC"], summary)

wilcox.test(mito_fraction ~ metadata_pseudobulk[, "ADNC"])

plot_df <- data.frame(mito_fraction = mito_fraction,
                      ADNC = metadata_pseudobulk[, "ADNC"])

gg <- ggplot2::ggplot(plot_df,
                      ggplot2::aes(x = ADNC, y = mito_fraction)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6) +
  ggplot2::labs(y = "Mitochondrial fraction of pseudobulk counts",
                title = "Ambient contamination is not balanced across groups") +
  ggplot2::theme_bw()
plot(gg)

dds_adj <- DESeq2::DESeqDataSetFromMatrix(
  countData = mat_pseudobulk,
  colData = metadata_pseudobulk,
  design = ~ sex + APOE4.status + mito_fraction + ADNC)

dds_adj <- DESeq2::DESeq(dds_adj)
DESeq2::resultsNames(dds_adj)

deseq2_res_adj <- DESeq2::results(dds_adj, name = "ADNC_AD_vs_NonAD")
summary(deseq2_res_adj)

deseq2_res_adj_ordered <- deseq2_res_adj[order(deseq2_res_adj$padj), ]
head(as.data.frame(deseq2_res_adj_ordered), 10)

sum(deseq2_res$padj < 0.05, na.rm = TRUE)
sum(deseq2_res_adj$padj < 0.05, na.rm = TRUE)

sum(grepl("^MT-", rownames(deseq2_res_adj_ordered)[1:20]))

# Why this had to be pseudobulk ------------------------------------------------

Seurat::Idents(seurat_obj) <- "ADNC"
adnc_binary <- rep("NonAD", ncol(seurat_obj))
adnc_binary[seurat_obj$ADNC %in% c("High", "Intermediate")] <- "AD"
seurat_obj$adnc_binary <- adnc_binary
Seurat::Idents(seurat_obj) <- "adnc_binary"

cell_level_res <- Seurat::FindMarkers(seurat_obj,
                                      ident.1 = "AD",
                                      ident.2 = "NonAD",
                                      test.use = "wilcox",
                                      logfc.threshold = 0,
                                      min.pct = 0.1)
head(cell_level_res, 5)

sum(cell_level_res$p_val_adj < 0.05, na.rm = TRUE)
sum(deseq2_res$padj < 0.05, na.rm = TRUE)

min(cell_level_res$p_val, na.rm = TRUE)
min(deseq2_res$pvalue, na.rm = TRUE)

shared_genes <- intersect(rownames(cell_level_res), rownames(deseq2_res))
length(shared_genes)

compare_df <- data.frame(
  gene = shared_genes,
  cell_level = -log10(cell_level_res[shared_genes, "p_val"] + 1e-300),
  pseudobulk = -log10(deseq2_res[shared_genes, "pvalue"]))

gg <- ggplot2::ggplot(compare_df,
                      ggplot2::aes(x = pseudobulk, y = cell_level)) +
  ggplot2::geom_point(alpha = 0.2, size = 0.6) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
  ggplot2::labs(x = "-log10 p, pseudobulk (donors as replicates)",
                y = "-log10 p, cell-level Wilcoxon (cells as replicates)",
                title = "The same genes, the same data, two definitions of n") +
  ggplot2::theme_bw()
plot(gg)

# How much of the variation is what you care about? ----------------------------

varpart_formula <- ~ ADNC + sex + APOE4.status + mito_fraction + n_cells

vsd_mat_top <- vsd_mat[order(apply(vsd_mat, 1, var),
                             decreasing = TRUE)[1:2000], ]

varpart_res <- variancePartition::fitExtractVarPartModel(
  vsd_mat_top, varpart_formula, metadata_pseudobulk)

head(variancePartition::sortCols(varpart_res))

sex_genes <- intersect(c("XIST", "UTY", "USP9Y", "RPS4Y1", "DDX3Y", "KDM5D"),
                       rownames(varpart_res))
varpart_res[sex_genes, ]

gg <- variancePartition::plotVarPart(variancePartition::sortCols(varpart_res)) +
  ggplot2::ggtitle("Variance explained per gene, 2000 most variable")
plot(gg)

round(100 * colMeans(varpart_res), 2)
round(100 * apply(varpart_res, 2, median), 2)

# What do the genes do? --------------------------------------------------------

teststat_vec <- deseq2_res[, "stat"]
names(teststat_vec) <- rownames(deseq2_res)
teststat_vec <- teststat_vec[!is.na(teststat_vec)]
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

length(teststat_vec)
head(teststat_vec)
tail(teststat_vec)

sum(grepl("^ENSG", names(teststat_vec)))

gse <- clusterProfiler::gseGO(teststat_vec,
                              ont = "BP",
                              keyType = "SYMBOL",
                              OrgDb = "org.Hs.eg.db",
                              minGSSize = 10,
                              maxGSSize = 500,
                              seed = TRUE,
                              verbose = FALSE)

gse_df <- as.data.frame(gse)
dim(gse_df)
head(gse_df[, c("Description", "setSize", "NES", "pvalue", "p.adjust")], 10)

gse_again <- clusterProfiler::gseGO(teststat_vec,
                                    ont = "BP",
                                    keyType = "SYMBOL",
                                    OrgDb = "org.Hs.eg.db",
                                    minGSSize = 10,
                                    maxGSSize = 500,
                                    seed = TRUE,
                                    verbose = FALSE)

identical(as.data.frame(gse_again)$pvalue, gse_df$pvalue)

gg <- enrichplot::dotplot(gse, showCategory = 15) +
  ggplot2::ggtitle("GO biological process, ranked by Wald statistic")
plot(gg)

gg <- enrichplot::dotplot(gse, showCategory = 10, split = ".sign") +
  ggplot2::facet_grid(. ~ .sign)
plot(gg)
