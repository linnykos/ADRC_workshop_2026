# Source of this tutorial: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
rm(list=ls())
# install.packages("SoupX")
library(SoupX)
library(Matrix)

# there's a few things you need to do:
# 1) make sure there's a folder called "filtered_feature_bc_matrix" and "raw_feature_bc_matrix"
# 2) both folders need the files "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz"
# 3) If there is an "analysis" folder, SoupX::load10X() will automatically grab both the clusters and the dimension reduction

# Load data and estimate soup profile
sc <- SoupX::load10X("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc_v2/")

library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
dd$clusters <- as.factor(paste0("Cluster_", dd$clusters))
mids = aggregate(cbind(tSNE1, tSNE2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = clusters), size = 0.2) + 
  geom_label(data = mids, aes(label = clusters)) + ggtitle("PBMC 1k Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

sc = SoupX::autoEstCont(sc)
out = SoupX::adjustCounts(sc)

cntSoggy = Matrix::rowSums(sc$toc > 0)
cntStrained = Matrix::rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(Matrix::rowSums(sc$toc > out)/Matrix::rowSums(sc$toc > 0)), n = 20)

plotChangeMap(sc, out, "IGKC")
