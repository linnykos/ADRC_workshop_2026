# Week 5: Bulk deconvolution with MuSiC
# Kevin Z. Lin, 2026-08-24
#
# This is the code from music-deconvolution_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data comes from the SeuratData package. Uncomment the InstallData()
# and install_github() lines the first time you run this; they download
# a few GB into your R library and only need to be run once.

# Setting up -------------------------------------------------------------------

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("TOAST")

# install.packages("devtools")
# devtools::install_github("xuranw/MuSiC")

# devtools::install_github("satijalab/seurat-data")

library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)

rm(list = ls())

# The data ---------------------------------------------------------------------

# SeuratData::InstallData("pbmcsca")

obj <- SeuratData::LoadData("pbmcsca")
obj

table(obj$CellType)

## Choosing which cell types to keep
select_ct <- c("B cell", "CD14+ monocyte", "CD4+ T cell",
               "Cytotoxic T cell", "Natural killer cell")

## Restricting to one technology
keep_vec <- obj$Method == "10x Chromium (v3)" & obj$CellType %in% select_ct
sum(keep_vec)

count_mat     <- SeuratObject::LayerData(obj, assay = "RNA",
                                         layer = "counts")[, keep_vec]
cell_type_vec <- obj$CellType[keep_vec]
dim(count_mat)
table(cell_type_vec)

# Splitting the cells into a reference half and a bulk half --------------------

set.seed(1)
ref_cell_vec  <- sample(colnames(count_mat), size = round(ncol(count_mat) / 2))
pool_cell_vec <- setdiff(colnames(count_mat), ref_cell_vec)

length(ref_cell_vec)
length(pool_cell_vec)
table(cell_type_vec[pool_cell_vec])

## The reference needs a notion of "subject"
set.seed(2)
donor_vec <- sample(paste0("donor", 1:4), size = length(ref_cell_vec),
                    replace = TRUE)

ref_meta_df <- data.frame(CellType  = cell_type_vec[ref_cell_vec],
                          sample_id = donor_vec,
                          row.names = ref_cell_vec)
table(ref_meta_df$sample_id, ref_meta_df$CellType)

## Packaging the reference
sce_ref <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = count_mat[, ref_cell_vec]),
  colData = ref_meta_df)
sce_ref

# Manufacturing the bulk samples -----------------------------------------------

n_sample <- 10

set.seed(10)
n_cell_mat <- matrix(sample(20:150, size = n_sample * length(select_ct),
                            replace = TRUE),
                     nrow = n_sample, ncol = length(select_ct),
                     dimnames = list(paste0("sample", 1:n_sample), select_ct))
n_cell_mat
rowSums(n_cell_mat)

pool_count_mat <- count_mat[, pool_cell_vec]
pool_ct_vec    <- cell_type_vec[pool_cell_vec]

set.seed(10)
bulk_mat <- sapply(rownames(n_cell_mat), function(sample_name) {
  drawn_vec <- unlist(lapply(select_ct, function(ct) {
    sample(pool_cell_vec[pool_ct_vec == ct],
           size = n_cell_mat[sample_name, ct], replace = TRUE)
  }))
  Matrix::rowSums(pool_count_mat[, drawn_vec])
})

dim(bulk_mat)
bulk_mat[1:5, 1:4]
colSums(bulk_mat)

prop_true_mat <- n_cell_mat / rowSums(n_cell_mat)
round(prop_true_mat, 3)

# Running MuSiC ----------------------------------------------------------------

music_res <- MuSiC::music_prop(bulk.mtx = bulk_mat, sc.sce = sce_ref,
                               clusters = "CellType", samples = "sample_id",
                               select.ct = select_ct)

names(music_res)

round(music_res$Est.prop.weighted, 3)

round(music_res$r.squared.full, 4)

# Scoring the estimates against the truth --------------------------------------

eval_mat <- MuSiC::Eval_multi(prop.real = prop_true_mat,
                              prop.est  = list(music_res$Est.prop.weighted,
                                               music_res$Est.prop.allgene),
                              method.name = c("MuSiC", "NNLS"))
eval_mat

plot_df <- data.frame(
  truth     = c(prop_true_mat, prop_true_mat),
  estimate  = c(music_res$Est.prop.weighted[, select_ct],
                music_res$Est.prop.allgene[, select_ct]),
  cell_type = factor(rep(rep(select_ct, each = n_sample), times = 2),
                     levels = select_ct),
  method    = factor(rep(c("MuSiC", "NNLS"),
                         each = n_sample * length(select_ct)),
                     levels = c("MuSiC", "NNLS")))

gg <- ggplot2::ggplot(plot_df, ggplot2::aes(x = truth, y = estimate)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "gray60",
                       linetype = "dashed") +
  ggplot2::geom_point(ggplot2::aes(colour = method), size = 2, alpha = 0.8) +
  ggplot2::facet_wrap(~ cell_type, nrow = 1) +
  ggplot2::coord_fixed(xlim = c(0, 0.55), ylim = c(0, 0.55)) +
  ggplot2::labs(x = "True proportion", y = "Estimated proportion",
                colour = "Method") +
  ggplot2::theme_minimal()
plot(gg)

err_music_mat <- abs(music_res$Est.prop.weighted[, select_ct] - prop_true_mat)
err_nnls_mat  <- abs(music_res$Est.prop.allgene[, select_ct]  - prop_true_mat)

round(rbind(MuSiC = colMeans(err_music_mat),
            NNLS  = colMeans(err_nnls_mat)), 4)

bar_df <- data.frame(
  proportion = c(prop_true_mat,
                 music_res$Est.prop.weighted[, select_ct],
                 music_res$Est.prop.allgene[, select_ct]),
  cell_type  = factor(rep(rep(select_ct, each = n_sample), times = 3),
                      levels = select_ct),
  sample     = factor(rep(rownames(prop_true_mat),
                          times = 3 * length(select_ct)),
                      levels = rownames(prop_true_mat)),
  source     = factor(rep(c("Truth", "MuSiC", "NNLS"),
                          each = n_sample * length(select_ct)),
                      levels = c("Truth", "MuSiC", "NNLS")))

gg <- ggplot2::ggplot(bar_df, ggplot2::aes(x = sample, y = proportion,
                                           fill = cell_type)) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::facet_wrap(~ source, ncol = 1) +
  ggplot2::labs(x = NULL, y = "Proportion", fill = "Cell type") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
plot(gg)

# What the weighting actually does ---------------------------------------------

weight_vec <- music_res$Weight.gene[, 1]
weight_vec <- weight_vec[!is.na(weight_vec)]
signif(quantile(weight_vec), 3)

names(sort(weight_vec))[1:10]

abundance_vec <- bulk_mat[names(weight_vec), 1] / sum(bulk_mat[, 1])
weight_df     <- data.frame(abundance = abundance_vec, weight = weight_vec)
weight_df     <- weight_df[weight_df$abundance > 0, ]
nrow(weight_df)
round(cor(log10(weight_df$abundance), log10(weight_df$weight)), 3)

gg <- ggplot2::ggplot(weight_df, ggplot2::aes(x = abundance, y = weight)) +
  ggplot2::geom_point(size = 0.4, alpha = 0.2) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "Gene's share of the bulk sample (log scale)",
                y = "MuSiC weight (log scale)") +
  ggplot2::theme_minimal()
plot(gg)

# Breaking it: a reference from the wrong assay --------------------------------

v2_vec      <- c("10x Chromium (v2)", "10x Chromium (v2) A",
                 "10x Chromium (v2) B")
keep_v2_vec <- obj$Method %in% v2_vec & obj$CellType %in% select_ct
sum(keep_v2_vec)
table(obj$Method[keep_v2_vec])

sce_ref_v2 <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = SeuratObject::LayerData(
                   obj, assay = "RNA", layer = "counts")[, keep_v2_vec]),
  colData = data.frame(CellType  = obj$CellType[keep_v2_vec],
                       sample_id = obj$Method[keep_v2_vec],
                       row.names = colnames(obj)[keep_v2_vec]))

music_v2_res <- MuSiC::music_prop(bulk.mtx = bulk_mat, sc.sce = sce_ref_v2,
                                  clusters = "CellType", samples = "sample_id",
                                  select.ct = select_ct, verbose = FALSE)

MuSiC::Eval_multi(prop.real = prop_true_mat,
                  prop.est  = list(music_res$Est.prop.weighted,
                                   music_v2_res$Est.prop.weighted),
                  method.name = c("MuSiC, v3 reference",
                                  "MuSiC, v2 reference"))

round(music_v2_res$Est.prop.weighted, 3)

round(colMeans(abs(music_v2_res$Est.prop.weighted[, select_ct]
                   - prop_true_mat)), 4)

tapply(obj$nCount_RNA, obj$Method, median)[c(v2_vec, "10x Chromium (v3)")]

round(music_v2_res$r.squared.full, 4)
