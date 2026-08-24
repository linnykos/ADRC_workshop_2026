# Week 4: Batch correction, part 2 — the server half
# Kevin Z. Lin, 2026-08-20
#
# This runs the slow part of crossdataset-integration_tutorial.Rmd / .html:
# finding anchors between SEA-AD and ROSMAP, correcting the expression matrix,
# and transferring labels from one study to the other. It does nothing else.
# See the .html for what each step is doing and why.
#
# Unlike the other scripts in this folder, this one is NOT meant to be read
# line by line in RStudio. It is submitted as a batch job:
#
#     sbatch crossdataset-integration_server.slurm
#
# It writes one file, crossdataset-integration_results.RData, which you then
# copy back to your own machine and load from the .Rmd. From your laptop:
#
#     scp you@server:/path/to/crossdataset-integration_results.RData .
#
# It downloads the two datasets itself (about 2 GB total) but it does NOT
# install anything — it checks that the packages are already there and stops
# immediately if they are not.
#
# Cost, measured: about 9 minutes end to end, of which FindIntegrationAnchors()
# is 4.4 min, IntegrateData() 0.3 min, the label transfer 0.2 min, and most of
# the rest is reading the 1.5 GB ROSMAP file. It writes a 376 MB .RData.
#
# That is small, and deliberately so -- see the tutorial's "Running the slow step
# somewhere else" for why a five-minute job is the right one to learn sbatch on.
# For contrast, the same script keeping all six ROSMAP brain regions (133,007
# query nuclei instead of 9,354) took 40 minutes and wrote 1.7 GB. Both calls
# scale with the number of query cells, and on a big enough query it is the
# memory rather than the time that bites. The script times the slow calls and
# prints the result, so after one run on YOUR data you will know what to ask for.

# Settings ---------------------------------------------------------------------

# Where to put the two downloaded datasets. tempdir() works, but on a cluster
# the node's /tmp can be small and 2 GB of downloads will overrun it. If the job
# dies during the download, point this at your scratch space instead, e.g.:
# data_dir <- "/gscratch/scrubbed/YOUR_USERNAME/adrc"
data_dir <- tempdir(check = TRUE)

# Where to write the result. The default is the directory you submitted from.
out_dir <- getwd()

results_file <- file.path(out_dir, "crossdataset-integration_results.RData")

# The three knobs the tutorial discusses. target_region is the ROSMAP brain
# region to keep: the shipped file holds all six that Sun et al. profiled, and
# pairing SEA-AD's middle temporal gyrus against any other one would make
# "study" and "brain region" the same variable. n_cells is how many nuclei to
# keep from SEA-AD, so that neither side dominates the anchor search; n_dims is
# how many dimensions the anchor search and the correction are allowed to see.
target_region <- "MidtemporalCortex"
n_cells <- 10000
n_dims <- 30

# Checking the packages are there ----------------------------------------------

# A batch job that waits an hour in the queue and then dies because one package
# is missing has cost you that hour. Check first, fail in the first second.
#
# This is the whole list, and it is shorter than the tutorial's: the .Rmd also
# needs RANN, ggplot2 and patchwork, and none of them are here because none of
# the plotting or the mixing score happens on the server.
required_packages <- c("Seurat", "SeuratObject")

missing_packages <- required_packages[!vapply(required_packages,
                                              requireNamespace,
                                              FUN.VALUE = logical(1),
                                              quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       ".\nInstall them on the server first, then resubmit:\n",
       '  install.packages(c("',
       paste(missing_packages, collapse = '", "'), '"))',
       call. = FALSE)
}

library(Seurat)
library(SeuratObject)

cat("Seurat", as.character(utils::packageVersion("Seurat")), "\n")
cat(R.version.string, "\n")
cat("RNGkind:", paste(RNGkind(), collapse = " / "), "\n")

# Getting the data -------------------------------------------------------------

# R's default download timeout is 60 seconds, which truncates a large file and
# reports it as a *warning* rather than an error, so the script would keep going
# and fail later at load(). The Dropbox links must end in dl=1; dl=0 serves an
# HTML preview page. mode = "wb" avoids silent corruption on Windows.
options(timeout = 3600)

seaad_url <- paste0("https://www.dropbox.com/scl/fi/cei02rwh3gx61uu9zzeql/",
                    "seaad_microglia.RData",
                    "?rlkey=04erfwn3w2lrnr80pefwfh8fu&dl=1")
seaad_file <- file.path(data_dir, "seaad_microglia.RData")
seaad_size <- 525672784

if (!file.exists(seaad_file) || file.size(seaad_file) != seaad_size) {
  download.file(seaad_url, destfile = seaad_file, mode = "wb")
}

rosmap_url <- paste0("https://www.dropbox.com/scl/fi/siecgiljoks48pkin42ts/",
                     "rosmap_microglia.RData",
                     "?rlkey=vwnmoeynhbsuxkiumz36hy301&dl=1")
rosmap_file <- file.path(data_dir, "rosmap_microglia.RData")
rosmap_size <- 1507998516

if (!file.exists(rosmap_file) || file.size(rosmap_file) != rosmap_size) {
  download.file(rosmap_url, destfile = rosmap_file, mode = "wb")
}

# Stop here rather than at load(), where the error would be about a corrupt
# connection and would not mention the download at all.
stopifnot(file.size(seaad_file) == seaad_size,
          file.size(rosmap_file) == rosmap_size)

load(seaad_file)
load(rosmap_file)

seaad <- seurat_obj
rm(list = c("seurat_obj")); gc(TRUE)

# Making the two objects comparable --------------------------------------------

# Everything in this section is a copy of the corresponding chunks in the .Rmd,
# in the same order, with the same seed. It has to be: the .Rmd runs this half
# on your laptop and then loads the object built here, so the two have to agree
# on which nuclei and which genes they are talking about. The .Rmd checks that
# they do, and that check is the only thing standing between you and a silent
# mismatch.
#
# The one intended difference: target_region, n_cells and n_dims are settings up
# top here and are written out as "MidtemporalCortex", 10000 and 1:30 in the
# .Rmd. Change one and you must change the other. The .Rmd's check will catch it
# if you forget -- all three travel back inside run_info for exactly that reason.

# The ROSMAP file ships all six brain regions. Cutting it down to one is the
# first thing that happens to it, before anything is renamed, subset by gene or
# normalized, so that every count printed below counts what we actually
# integrate. Assert the region exists rather than discovering an empty object
# forty minutes later.
stopifnot(target_region %in% rosmap$brainRegion)

rosmap_cells <- Seurat::Cells(rosmap)[rosmap$brainRegion == target_region]
rosmap <- subset(rosmap, cells = rosmap_cells)

cat("ROSMAP", target_region, "nuclei:", ncol(rosmap), "\n")
cat("ROSMAP donors:", length(unique(rosmap$subject)), "\n")

seaad <- Seurat::RenameCells(seaad, add.cell.id = "seaad")
rosmap <- Seurat::RenameCells(rosmap, add.cell.id = "rosmap")

seaad$dataset <- "SEA-AD"
rosmap$dataset <- "ROSMAP"

common_genes <- intersect(rownames(seaad), rownames(rosmap))
cat("common genes:", length(common_genes), "\n")

seaad <- subset(seaad, features = common_genes)
rosmap <- subset(rosmap, features = common_genes)

set.seed(10)
seaad_cells <- sample(Seurat::Cells(seaad), size = n_cells)
seaad <- subset(seaad, cells = seaad_cells)

cat("cells:", ncol(seaad), "SEA-AD,", ncol(rosmap), "ROSMAP\n")

seaad <- Seurat::NormalizeData(seaad, verbose = FALSE)
rosmap <- Seurat::NormalizeData(rosmap, verbose = FALSE)

seaad <- Seurat::FindVariableFeatures(seaad, selection.method = "vst",
                                      nfeatures = 2000, verbose = FALSE)
rosmap <- Seurat::FindVariableFeatures(rosmap, selection.method = "vst",
                                       nfeatures = 2000, verbose = FALSE)

# Allow up to 16 GB (16 * 1024^3 bytes) for exported globals
options(future.globals.maxSize = 16 * 1024^3)
integration_features <- Seurat::SelectIntegrationFeatures(
  object.list = list(seaad, rosmap),
  nfeatures = 2000,
  verbose = FALSE)

cat("integration features:", length(integration_features), "\n")

# Integrating ------------------------------------------------------------------

# The two calls this whole script exists for. Timing each one, because the .out
# log is the only place you will find out where the time actually went.
time_start <- Sys.time()

integration_anchors <- Seurat::FindIntegrationAnchors(
  object.list = list(seaad, rosmap),
  anchor.features = integration_features,
  reduction = "cca",
  dims = 1:n_dims,
  verbose = FALSE)

time_anchors <- Sys.time()
cat("FindIntegrationAnchors:",
    round(difftime(time_anchors, time_start, units = "mins"), 1), "min\n")

integrated <- Seurat::IntegrateData(anchorset = integration_anchors,
                                    new.assay.name = "integrated",
                                    dims = 1:n_dims,
                                    verbose = FALSE)

time_integrate <- Sys.time()
cat("IntegrateData:",
    round(difftime(time_integrate, time_anchors, units = "mins"), 1), "min\n")

# Note what is NOT done here: no JoinLayers(), no ScaleData(), no RunPCA(), no
# RunUMAP(). Those are seconds of work and they belong in the .Rmd, where you
# can see their output and change your mind about them without resubmitting a
# job. Only put the expensive, settled steps on the server.

# Transferring labels ----------------------------------------------------------

# The alternative route, which never edits an expression value. Same anchor
# machinery underneath, so it is slow for the same reason and belongs here.
transfer_anchors <- Seurat::FindTransferAnchors(reference = seaad,
                                                query = rosmap,
                                                dims = 1:n_dims,
                                                verbose = FALSE)

predictions_df <- Seurat::TransferData(anchorset = transfer_anchors,
                                       refdata = seaad$Supertype,
                                       dims = 1:n_dims,
                                       verbose = FALSE)

cat("TransferData:",
    round(difftime(Sys.time(), time_integrate, units = "mins"), 1), "min\n")

# Saving -----------------------------------------------------------------------

# run_info is what makes the reproducibility check in the .Rmd mean something.
# The .Rmd rebuilds seaad and rosmap on your laptop and then loads an object
# built here; that only lines up if both machines drew the same 10,000 nuclei,
# which in turn depends on the RNG being the one R has used by default since
# 3.6.0. Record it rather than assume it.
run_info <- list(r_version = R.version.string,
                 seurat_version = as.character(utils::packageVersion("Seurat")),
                 rng_kind = RNGkind(),
                 target_region = target_region,
                 n_cells = n_cells,
                 n_dims = n_dims,
                 run_at = Sys.time())

save(integrated, predictions_df, seaad_cells, integration_features, run_info,
     file = results_file)

cat("wrote", results_file, "\n")
cat("size:", round(file.size(results_file) / 1e6), "MB\n")
cat("copy it back with:\n")
cat("  scp", paste0("USER@SERVER:", results_file), ".\n")
