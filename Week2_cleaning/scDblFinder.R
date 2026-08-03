rm(list=ls())
library(Seurat)
library(scDblFinder)

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

clustering_df <- read.csv(file.path(tmp_dir, "analysis", "clustering",
                                    "graphclust", "clusters.csv"), row.names = 1)
umap_df <- read.csv(file.path(tmp_dir, "analysis", "umap",
                                    "2_components", "projection.csv"), row.names = 1)
umap_df <- as.matrix(umap_df)

pbmc.data <- Seurat::Read10X(data.dir = paste0(data_dir, "/filtered_feature_bc_matrix"))

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data[["Gene Expression"]], 
                           project = "pbmc1k", 
                           meta.data = clustering_df)
pbmc
pbmc[["umap"]] <- Seurat::CreateDimReducObject(umap_df, key = "RNA")

# https://satijalab.org/seurat/reference/as.singlecellexperiment
Seurat::DefaultAssay(pbmc) <- "RNA"

# Only needed if you have Seurat version 5
rna_counts <- Seurat::GetAssayData(pbmc[["RNA"]], layer = "counts")
pbmc[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_counts)

sce <- Seurat::as.SingleCellExperiment(pbmc)

# This procedure is a bit "better" if you have cell cluster information
## This code below though is simply for pedagogical demonstration
# https://bioconductor.org/books/3.15/OSCA.advanced/doublet-detection.html
set.seed(10)
sce.dbl <- scDblFinder::scDblFinder(sce, clusters=sce$Cluster) # takes about 1 minute

# Take a look at the doublet scores
head(sce.dbl@colData$scDblFinder.score)
quantile(sce.dbl@colData$scDblFinder.score)
hist(sce.dbl@colData$scDblFinder.score)

# scDblFinder also enumerates the cells with "too high" of a doublet score
table(sce.dbl$scDblFinder.class)

pbmc$scDblFinder.score <- sce.dbl@colData$scDblFinder.score
pbmc$scDblFinder.class <- sce.dbl$scDblFinder.class

Seurat::DimPlot(pbmc, group.by = "scDblFinder.class")

pbmc <- subset(pbmc, scDblFinder.class == "singlet")
pbmc
