#!/usr/bin/env Rscript
# Stage 3 prep: export the QC-passing barcode list per sample (from seurat_qc/<SRR>.rds)
# so scrublet_doublets.py can restrict its input matrix to the same cells that
# survived Stage 2 filtering.

suppressPackageStartupMessages(library(Seurat))

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)

runs <- readLines("runs.txt")

for (srr in runs) {
  rds_path <- file.path("seurat_qc", paste0(srr, ".rds"))
  if (!file.exists(rds_path)) {
    warning("Missing ", rds_path, " -- skipping")
    next
  }
  obj <- readRDS(rds_path)
  writeLines(colnames(obj), file.path("seurat_qc", paste0(srr, "_barcodes.txt")))
  message(srr, ": ", ncol(obj), " QC-passing barcodes written")
}
