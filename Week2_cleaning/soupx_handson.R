# Source of this tutorial: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
rm(list=ls())
# install.packages("SoupX")
library(SoupX)
library(Matrix)
library(ggplot2)

# there's a few things you need to do:
# 1) make sure there's a folder called "filtered_feature_bc_matrix" and "raw_feature_bc_matrix"
# 2) both folders need the files "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz"
# 3) If there is an "analysis" folder, SoupX::load10X() will automatically grab both the clusters and the dimension reduction

# Load data and estimate soup profile
tmpDir = tempdir(check = TRUE)
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "tod.tar.gz"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "toc.tar.gz"))
untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)

toc = Seurat::Read10X(file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38"))
tod = Seurat::Read10X(file.path(tmpDir, "raw_gene_bc_matrices", "GRCh38"))
sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc = estimateSoup(sc)

# You could need to write different lines to set the cluster and 2-component dimension reduction (i.e. UMAP) here with slightly different code
data(PBMC_metaData)
rownames(PBMC_metaData) <- paste0(rownames(PBMC_metaData), "-1")
sc = setClusters(sc, setNames(PBMC_metaData$Annotation, rownames(PBMC_metaData)))
sc = setDR(sc, PBMC_metaData[colnames(sc$toc), c("RD1", "RD2")])

dd = sc$metaData[colnames(sc$toc),]
mids = aggregate(cbind(RD1, RD2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = clusters), size = 0.2) + 
  geom_label(data = mids, aes(label = clusters)) + ggtitle("PBMC 4k Clusters") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

gg = plotMarkerMap(sc, "IGKC")
plot(gg)

##########

marker_genes <- SoupX::quickMarkers(sc$toc, sc$metaData[,"clusters"])
cluster_genes <- split(marker_genes$gene, marker_genes$cluster)
cluster_genes

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = cluster_genes)
for(i in 1:length(cluster_genes)){
  gg = plotMarkerMap(sc, geneSet = cluster_genes[[i]], useToEst = useToEst[,i])
  gg <- gg + ggplot2::ggtitle(paste0("Cluster ", i))
  plot(gg)
}

sc = calculateContaminationFraction(sc, cluster_genes, useToEst = useToEst)
head(sc$metaData)
out = SoupX::adjustCounts(sc)

cntSoggy = Matrix::rowSums(sc$toc > 0)
cntStrained = Matrix::rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(Matrix::rowSums(sc$toc - out)), n = 20)

plotChangeMap(sc, out, "RPL10")
plotChangeMap(sc, out, "LYN")
