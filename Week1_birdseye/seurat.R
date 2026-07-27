# https://satijalab.org/seurat/articles/seurat5_essential_commands 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial 

rm(list=ls())
library(Seurat)
pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k")
pbmc

# look at cell-level metadata
pbmc@meta.data

# extract the count matrix. It's a sparse matrix, 
## which you can do matrix operations, as if any other matrix in R
count_mat <- SeuratObject::LayerData(pbmc,
                                     layer = "counts",
                                     assay = "RNA")
count_mat[42:44,14:16]
head(colSums(count_mat)) # the total number of counts per cell
head(pbmc$nCount_RNA) # notice the previous line gets you the same numbers as this line

# extract the gene names
SeuratObject::Features(pbmc)[1:5]

# extract the cell names
Seurat::Cells(pbmc)[1:5]

# Quality control
# 1) Mitochondrial percentage
# First see: how many mitochondrial genes are there?
gene_vec <- SeuratObject::Features(pbmc)
gene_vec[grep("^MT-", gene_vec)]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, 
                                                     pattern = "^MT-")

# 2) Ribosomal percentages
gene_vec[grep("^RPS", gene_vec)]
pbmc[["percent.rb"]] <- Seurat::PercentageFeatureSet(pbmc, 
                                                     pattern = "^RPS")

# 3) Cell Cycling
# see more: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
pbmc <- Seurat::NormalizeData(pbmc) #Seurat::CellCycleScoring requires the "data" slot filled
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- Seurat::CellCycleScoring(pbmc, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes)

head(pbmc@meta.data)

Seurat::VlnPlot(pbmc, features = c("nFeature_RNA",  "nCount_RNA", "percent.mt", "percent.rb"), 
                ncol = 4)

pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

