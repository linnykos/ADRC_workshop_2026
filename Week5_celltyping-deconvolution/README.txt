Slides for Week 5: TBD

music-deconvolution_tutorial.Rmd / .html / .R
  Bulk deconvolution with MuSiC. Builds ten synthetic bulk samples by summing
  counts over known mixtures of pbmcsca cells, deconvolves them against a
  single-cell reference held out from the same cells, and scores the estimates
  against the mixtures we chose. Closes by swapping in a reference from a
  different 10x chemistry and watching the estimates collapse.

  Runs in about 15 seconds on a laptop. The only data it needs is pbmcsca from
  SeuratData, which Week 3 already installed.
