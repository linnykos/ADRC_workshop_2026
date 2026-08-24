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
      no shared donors, matched brain region. Subsetting ROSMAP to one region
      and reconciling barcodes/genes/size, a nearest-neighbour batch-mixing
      score measured before and after, FindIntegrationAnchors() +
      IntegrateData(), what is actually inside the "integrated" assay (no
      counts, ~2,000 genes, negative values, zeros filled in), WHICH of the two
      studies actually got corrected, why you must not test genes on it, and
      label transfer as the alternative that never edits expression.
      Built from Kevin's Writeup30b exploration scripts.

      Measured on the render: uncorrected mixing median 0 for both studies;
      after correction SEA-AD 0.43 (null 0.48) but ROSMAP only 0.07 (null 0.52).
      And IntegrateData() left all 10,000 SEA-AD values BYTE-IDENTICAL while
      moving ROSMAP by up to 2.69 -- every negative value in the corrected
      matrix is ROSMAP's. Nothing in the call says which study is the reference;
      it falls out of the order of list(seaad, rosmap). That asymmetry is now a
      section of its own.

      THIS ONE IS SPLIT ACROSS TWO MACHINES, and the honest reason is
      pedagogical rather than computational. On this data the three anchor calls
      take about five minutes (FindIntegrationAnchors 4.4 min, IntegrateData
      0.3, the label transfer 0.2), so they are shown in the tutorial but marked
      eval=FALSE and run on a cluster instead -- because a five-minute job is
      the right size on which to get sbatch wrong three times. The same script
      keeping all six ROSMAP regions took 40 minutes and wrote 1.7 GB, and that
      contrast is in the tutorial. Running R as a batch job is a teaching goal
      of the tutorial in its own right. The two extra files:

  crossdataset-integration_server.R
      Self-contained batch script: downloads both datasets, repeats the data
      preparation, runs FindIntegrationAnchors() + IntegrateData() and the
      label transfer, and saves crossdataset-integration_results.RData. It
      does not install packages -- it checks they are present and stops with
      the install line if not. Nothing downstream (no scaling, PCA, UMAP or
      plots); those stay in the .Rmd where you can see their output.
      Not meant to be read line by line in RStudio, unlike everything else here.

  crossdataset-integration_server.slurm
      The SLURM job description for the above. You MUST fill in --account and
      --partition for your own cluster; the rest has working defaults
      (4 cores, 64G, 2 hours -- generous for the job as shipped; the
      six-region variant is what wanted 128G). Submit with:
          sbatch crossdataset-integration_server.slurm
      then copy the result back to sit next to the .Rmd:
          scp you@server:/path/to/crossdataset-integration_results.RData .
      The .Rmd stops with an instruction if it cannot find that file, and
      checks that the object inside it describes the same nuclei your own
      machine selected before it uses it.

DATA: the first two tutorials download seaad_microglia.RData (525 MB)
themselves into a temporary folder, so they run as-is with nothing to set up.
The third downloads that file plus rosmap_microglia.RData. The first run of
each spends about five minutes on the download. The SEA-AD file is a microglia/PVM
subset of SEA-AD middle temporal gyrus, downsampled to 40,000 nuclei, 36,412
genes, 89 donors. It arrives with counts and log-normalized data, the authors'
scVI integration, a UMAP, and 43 metadata columns including ADNC, Braak stage,
CERAD score, APOE4 status and donor ID. If you already have the file, point the
data_dir line near the top at the folder holding it and the download is skipped.

The ROSMAP file is the Sun et al., Cell 2023 microglia -- the WHOLE set, all six
brain regions, 133,007 nuclei from 425 individuals, 16,219 genes. It arrives
already processed: counts, data AND scale.data, its own pca/umap, and a
seurat_clusters column. Metadata columns are State (MG0-MG12, twelve of the
thirteen present), subject, brainRegion, batch, ADdiag3types, age_death, msex,
pmi, percent.mt, percent.rp.

  CAUTION: those are NOT the names an earlier draft of the tutorial assumed
  (state / donor_id / sequencing_batch / ad_diagnosis), and the object inside
  is named `rosmap`, not `rosmap_obj`. See private/data/rosmap/README.txt.

The tutorial and the server script both subset it to
brainRegion == "MidtemporalCortex" -- 9,354 nuclei, 48 individuals, a single
batch value -- as their first operation, and each records target_region so the
two halves can check they agree. The region was matched to SEA-AD's on purpose:
pairing SEA-AD MTG against ROSMAP's much larger PFC subset would have made
"study" and "brain region" the same variable, and keeping all six regions would
have made region vary INSIDE one of the two batches while leaving ROSMAP
outnumbering SEA-AD 13 to 1.

  Note also that ROSMAP's `batch` column is not a sequencing batch for most
  regions: five of the six carry the region's own name, and only PFC is
  resolved into fourteen SM_* runs.

  Both download links were checked again on 2026-08-24 by an actual knit and are
  live: the SEA-AD file returns 525,672,784 bytes and the ROSMAP file
  1,507,998,516 bytes, matching the sizes the scripts assert. As of 2026-08-24
  the third tutorial HAS been knitted and every number in its prose is read off
  that render.

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