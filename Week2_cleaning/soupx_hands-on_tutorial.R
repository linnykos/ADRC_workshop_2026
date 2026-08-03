# Week 2: Removing ambient RNA with SoupX (the manual route)
# Kevin Z. Lin, 2026-08-03
#
# This is the code from soupx_hands-on_tutorial.Rmd / .html, with nothing else.
# See the .html for what each step is doing and why.
#
# The data downloads itself into a temporary folder, so this script should
# run as-is. Nothing is written to your machine outside that folder.

# Setting up -------------------------------------------------------------------

library(SoupX)
library(Matrix)
library(ggplot2)

rm(list = ls())

# Getting the data -------------------------------------------------------------

tmp_dir <- tempdir(check = TRUE)
base_url <- "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k"

for (component in c("raw_gene_bc_matrices",
                    "filtered_gene_bc_matrices",
                    "analysis")) {
  destination <- file.path(tmp_dir, paste0(component, ".tar.gz"))
  download.file(paste0(base_url, "_", component, ".tar.gz"),
                destfile = destination)
  untar(destination, exdir = tmp_dir)
}

# Building the SoupChannel by hand ---------------------------------------------

toc <- Seurat::Read10X(file.path(tmp_dir, "filtered_gene_bc_matrices", "GRCh38"))
tod <- Seurat::Read10X(file.path(tmp_dir, "raw_gene_bc_matrices", "GRCh38"))

dim(toc)
dim(tod)

sc <- SoupX::SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- SoupX::estimateSoup(sc)
sc

head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)

# Attaching clusters and dimension reduction after the fact --------------------

clustering_df <- read.csv(file.path(tmp_dir, "analysis", "clustering",
                                    "graphclust", "clusters.csv"))
head(clustering_df)
table(clustering_df$Cluster)

sc <- SoupX::setClusters(sc, setNames(clustering_df$Cluster,
                                      clustering_df$Barcode))

tsne_df <- read.csv(file.path(tmp_dir, "analysis", "tsne",
                              "2_components", "projection.csv"))
head(tsne_df)

rownames(tsne_df) <- tsne_df$Barcode
sc <- SoupX::setDR(sc, tsne_df[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

head(sc$metaData)

# A first look at the data -----------------------------------------------------

dd <- sc$metaData[colnames(sc$toc), ]
dd$clusters <- as.factor(paste0("Cluster_", dd$clusters))

mids <- aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)

gg <- ggplot2::ggplot(dd, ggplot2::aes(TSNE.1, TSNE.2)) +
  ggplot2::geom_point(ggplot2::aes(colour = clusters), size = 0.2) +
  ggplot2::geom_label(data = mids, ggplot2::aes(label = clusters)) +
  ggplot2::ggtitle("PBMC 4k, CellRanger graph-based clusters") +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1)))
plot(gg)

# Choosing genes to estimate the contamination with ----------------------------

## Why the most abundant soup genes are the wrong choice
gg <- SoupX::plotMarkerMap(sc, "RPL10")
plot(gg)

## Letting the data suggest candidates
SoupX::plotMarkerDistribution(sc)

marker_genes <- SoupX::quickMarkers(sc$toc, sc$metaData$clusters)
head(marker_genes[, c("gene", "cluster", "tfidf", "qval")])

cluster_genes <- split(marker_genes$gene, marker_genes$cluster)
lengths(cluster_genes)

# Choosing which cells to estimate from ----------------------------------------

useToEst <- SoupX::estimateNonExpressingCells(sc, nonExpressedGeneList = cluster_genes)
dim(useToEst)

SoupX::plotMarkerMap(sc, geneSet = cluster_genes[[1]], useToEst = useToEst[, 1])
SoupX::plotMarkerMap(sc, geneSet = cluster_genes[[5]], useToEst = useToEst[, 5])

# Calculating the contamination fraction ---------------------------------------

sc <- SoupX::calculateContaminationFraction(sc, cluster_genes, useToEst = useToEst)
head(sc$metaData)

# Correcting the counts --------------------------------------------------------

out <- SoupX::adjustCounts(sc)

# What actually changed --------------------------------------------------------

cnt_soggy <- Matrix::rowSums(sc$toc > 0)
cnt_strained <- Matrix::rowSums(out > 0)
most_zeroed <- tail(sort((cnt_soggy - cnt_strained) / cnt_soggy), n = 10)
most_zeroed

tail(sort(Matrix::rowSums(sc$toc > out) / Matrix::rowSums(sc$toc > 0)), n = 20)

## The correction, seen directly
SoupX::plotChangeMap(sc, out, "LYZ")

SoupX::plotChangeMap(sc, out, "IGKC")

SoupX::plotChangeMap(sc, out, "RPL10")

# Using the cleaned matrix downstream ------------------------------------------

pbmc <- Seurat::CreateSeuratObject(counts = out)
