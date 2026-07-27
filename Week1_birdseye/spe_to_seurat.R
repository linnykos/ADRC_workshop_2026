#' Convert a SpatialExperiment to a Seurat object, expression only
#'
#' Ensembl rownames are always replaced by rowData(spe)$gene_name, because the
#' downstream workflow needs symbols: `PercentageFeatureSet(pattern = "^MT-")`
#' returns 0% for every spot on Ensembl IDs -- silently, with no error -- and
#' `CellCycleScoring()` matches nothing.
#'
#' @param spe          A SpatialExperiment (or SingleCellExperiment) object.
#' @param assay_name   Assay in `spe` to use as the Seurat "counts" layer.
#'                     Should be raw counts.
#' @param meta_cols    colData columns to carry over: "all" (default), "none",
#'                     or a character vector of column names.
#' @param sample_ids   Optional character vector; keep only these `colData(spe)$sample_id`
#'                     values. Handy for cutting 38k spots down to a laptop-sized subset.
#' @param project      Seurat project name.
#' @param min_cells,min_features Passed to `CreateSeuratObject()`. Defaults of 0 keep
#'                     everything, so filtering stays visible in the analysis script.
#'
#' @return A Seurat object with a single "RNA" assay and no spatial information.
spe_to_seurat <- function(spe,
                          assay_name = "counts",
                          meta_cols = "all",
                          sample_ids = NULL,
                          project = "spatial_expression",
                          min_cells = 0,
                          min_features = 0) {
  stopifnot(assay_name %in% SummarizedExperiment::assayNames(spe))

  ## Optional subsetting by sample, done before pulling the matrix out
  if (!is.null(sample_ids)) {
    keep_idx <- which(SummarizedExperiment::colData(spe)$sample_id %in% sample_ids)
    if (length(keep_idx) == 0) stop("No spots matched the requested sample_ids.")
    spe <- spe[, keep_idx]
  }

  col_data <- SummarizedExperiment::colData(spe)
  gene_vec <- as.character(SummarizedExperiment::rowData(spe)$gene_name)

  ## The two things that must hold for the renaming below to be safe
  stopifnot(!anyNA(gene_vec), all(gene_vec != ""))
  stopifnot("key" %in% colnames(col_data), anyDuplicated(col_data$key) == 0)

  ## The counts matrix, renamed: Ensembl IDs -> gene symbols on the rows, and
  ## barcode_sampleid on the columns. Raw Visium barcodes repeat across samples,
  ## so colnames(spe) cannot be used as Seurat cell names -- colData$key can.
  count_mat <- SummarizedExperiment::assay(spe, assay_name)
  rownames(count_mat) <- make.unique(gene_vec)   # 8 symbols are duplicated
  colnames(count_mat) <- as.character(col_data$key)

  ## Per-spot metadata. Everything here is non-spatial by construction --
  ## spatialCoords() and imgData() live in separate slots and are never touched.
  meta_df <- NULL
  if (!identical(meta_cols, "none")) {
    keep_cols <- if (identical(meta_cols, "all")) colnames(col_data) else meta_cols
    stopifnot(all(keep_cols %in% colnames(col_data)))
    meta_df <- as.data.frame(col_data[, keep_cols, drop = FALSE])
    rownames(meta_df) <- colnames(count_mat)
  }

  Seurat::CreateSeuratObject(counts = count_mat,
                             project = project,
                             meta.data = meta_df,
                             min.cells = min_cells,
                             min.features = min_features)
}
