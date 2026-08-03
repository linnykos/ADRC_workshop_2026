# Week 2: Removing ambient RNA with SoupX (the automated route)
# Kevin Z. Lin, 2026-08-03
#
# This is the code from soupx_tutorial.rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data downloads itself into a temporary folder, so this script should
# run as-is. Nothing is written to your machine outside that folder.

# Setting up -------------------------------------------------------------------

# install.packages("SoupX")

library(SoupX)
library(Matrix)
library(ggplot2)

rm(list = ls())

# Getting the data -------------------------------------------------------------

tmp_dir <- tempdir(check = TRUE)
base_url <- "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k"

for (component in c("raw_feature_bc_matrix",
                    "filtered_feature_bc_matrix",
                    "analysis")) {
  destination <- file.path(tmp_dir, paste0(component, ".tar.gz"))
  download.file(paste0(base_url, "_", component, ".tar.gz"),
                destfile = destination)
  untar(destination, exdir = tmp_dir)
}

data_dir <- tmp_dir

# If you would rather use your OWN CellRanger output, comment out the
# download loop above and point data_dir at the folder holding
# raw_feature_bc_matrix/, filtered_feature_bc_matrix/ and analysis/
# (usually the 'outs' directory of a CellRanger run), e.g.:
# data_dir <- "/path/to/your/cellranger/outs"

# Loading the data -------------------------------------------------------------

sc <- SoupX::load10X(data_dir)
sc

dim(sc$toc)
head(sc$metaData)

# Profiling the soup -----------------------------------------------------------

head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)

# A first look at the data -----------------------------------------------------

dd <- sc$metaData[colnames(sc$toc), ]
dd$clusters <- as.factor(paste0("Cluster_", dd$clusters))

mids <- aggregate(cbind(tSNE1, tSNE2) ~ clusters, data = dd, FUN = mean)

gg <- ggplot2::ggplot(dd, ggplot2::aes(tSNE1, tSNE2)) +
  ggplot2::geom_point(ggplot2::aes(colour = clusters), size = 0.2) +
  ggplot2::geom_label(data = mids, ggplot2::aes(label = clusters)) +
  ggplot2::ggtitle("PBMC 1k annotation") +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1)))
plot(gg)

## Why we suspect contamination in the first place
gg <- SoupX::plotMarkerMap(sc, "LYZ")
plot(gg)

# Estimating the contamination fraction ----------------------------------------

sc <- SoupX::autoEstCont(sc)

head(sc$metaData)

# Correcting the counts --------------------------------------------------------

out <- SoupX::adjustCounts(sc)

# What actually changed --------------------------------------------------------

cnt_soggy <- Matrix::rowSums(sc$toc > 0)
cnt_strained <- Matrix::rowSums(out > 0)
most_zeroed <- tail(sort((cnt_soggy - cnt_strained) / cnt_soggy), n = 10)
most_zeroed

tail(sort(Matrix::rowSums(sc$toc > out) / Matrix::rowSums(sc$toc > 0)), n = 20)

## Back to the lysozyme-producing T cells
SoupX::plotChangeMap(sc, out, "LYZ")

# Using the cleaned matrix downstream ------------------------------------------

pbmc <- Seurat::CreateSeuratObject(counts = out)
