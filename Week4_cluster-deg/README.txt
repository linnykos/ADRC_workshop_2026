Slides for Week 4: https://docs.google.com/presentation/d/1eDU3AZ3-vlEQ5aPBo1CCvck1-HNxUMqfB0eMDW42lPA/edit?usp=sharing

Week 4: Differential expression and cell-type labeling (plus batch correction, part 2)

Three tutorials. The first two run on the same dataset; the third adds a second
one alongside it.

  clustering-markers_tutorial.Rmd / .html / .R
      Clustering 40,000 microglial nuclei from 89 donors. Uncorrected PCA vs the
      SEA-AD scVI embedding, donor purity, cluster concordance, bootstrap
      stability, choosing a resolution, marker genes, curated marker panels, and
      a check against SEA-AD's own supertype labels.
      Built from the 2025 script clustering-markers.R.

  pseudobulk-deg_tutorial.Rmd / .html / .R
      Differential expression between donors. Pseudobulk aggregation, gene and
      sample filtering, design-matrix checks, DESeq2, apeglm shrinkage, why the
      top of the gene list is probably ambient RNA, a side-by-side against a
      cell-level test, variancePartition, and GSEA.
      Built from the 2025 scripts deseq2.R and gsea.R, merged.

  crossdataset-integration_tutorial.Rmd / .html / .R
      Batch correction, part 2. Week 3 produced an integrated *reduction* and
      stopped; this one produces a corrected *expression matrix* and then
      explains what it is and is not good for. SEA-AD middle temporal gyrus
      against ROSMAP middle temporal cortex (Sun et al. 2023) -- two studies,
      no shared donors, matched brain region. Reconciling barcodes/genes/size,
      a nearest-neighbour batch-mixing score measured before and after,
      FindIntegrationAnchors() + IntegrateData(), what is actually inside the
      "integrated" assay (no counts, ~2,000 genes, negative values, zeros
      filled in), why you must not test genes on it, and label transfer as the
      alternative that never edits expression.
      Built from Kevin's Writeup30b exploration scripts.

DATA: the first two tutorials download seaad_microglia.RData (525 MB)
themselves into a temporary folder, so they run as-is with nothing to set up.
The third downloads that file plus rosmap_microglia.RData. The first run of
each spends about five minutes on the download. The SEA-AD file is a microglia/PVM
subset of SEA-AD middle temporal gyrus, downsampled to 40,000 nuclei, 36,412
genes, 89 donors. It arrives with counts and log-normalized data, the authors'
scVI integration, a UMAP, and 43 metadata columns including ADNC, Braak stage,
CERAD score, APOE4 status and donor ID. If you already have the file, point the
data_dir line near the top at the folder holding it and the download is skipped.

The ROSMAP file is the middle temporal cortex microglia from Sun et al., Cell
2023 -- around 10,000 nuclei from 48 individuals, counts only, carrying the
published MG0-MG12 state labels plus donor, diagnosis, age, sex and PMI. The
region was matched to SEA-AD's on purpose: pairing SEA-AD MTG against ROSMAP's
much larger PFC subset would have made "study" and "brain region" the same
variable.

  ** NOT YET LIVE ** crossdataset-integration_tutorial.Rmd currently carries a
  PLACEHOLDER Dropbox link and a zero byte-count for rosmap_microglia.RData,
  and cannot knit until they are filled in. Build the file with
  private/code/rosmap/prep_rosmap_microglia_claude.R, upload it, and paste the
  dl=1 link and the printed byte count into the download_rosmap chunk (then
  regenerate the .R).

NOTE ON THE DOWNLOAD: R's default download.file() timeout is 60 seconds, which
is nowhere near enough for 525 MB. The tutorials set options(timeout = 3600) and
then verify the downloaded byte count. Without that, the transfer dies partway
through and R reports it as a WARNING rather than an error -- so a script can
sail past a truncated file and fail later at load() with a confusing message
about corruption.

NOTE ON KNITTING: all three documents set cache = TRUE plus autodep = TRUE. The first
knit of the clustering tutorial takes roughly 20 minutes (plus the download);
later knits take a couple of minutes. The _cache/ folders are large (GB scale)
and are gitignored. If a knit is interrupted, delete the _cache/ folder before 
re-running or you will get "error reading from connection". Delete _cache/ and 
knit once from clean before distributing, and check the HTML for empty images: 
pandoc's --embed-resources step has been seen to silently emit blank figures 
under disk pressure while exiting 0.