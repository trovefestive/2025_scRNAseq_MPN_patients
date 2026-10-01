#!/usr/bin/env Rscript
# Stage 3b: merge scrublet_doublets.py's per-cell calls back into each Stage-2
# Seurat object, drop predicted doublets (pre-merge, per plan), and stamp a
# qc_flag column so flagged low-quality samples are traceable through every
# downstream QC/cluster-composition plot.

suppressPackageStartupMessages(library(Seurat))

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)

runs <- readLines("runs.txt")

# From Stage 1/2 review (see notes/mpn-scrnaseq-project-status.md):
# SRR27748964 kept per user decision (2026-08-26) but flagged -- only 884
# post-QC-filter cells (34.8% retained) from a low median-genes/cell (160) library.
qc_flags <- list(
  SRR27748964 = "low_cell_recovery_884cells"
)

summary_rows <- list()

for (srr in runs) {
  rds_path <- file.path("seurat_qc", paste0(srr, ".rds"))
  scrub_path <- file.path("seurat_qc", paste0(srr, "_scrublet.csv"))
  if (!file.exists(rds_path) || !file.exists(scrub_path)) {
    warning("Missing input for ", srr, " -- skipping")
    next
  }

  obj <- readRDS(rds_path)
  scrub <- read.csv(scrub_path, stringsAsFactors = FALSE)
  rownames(scrub) <- scrub$barcode
  scrub <- scrub[colnames(obj), ]  # align order to Seurat object

  obj$doublet_score <- scrub$doublet_score
  # scrublet writes Python-style True/False, which read.csv keeps as character
  # (it only auto-coerces TRUE/FALSE) -- coerce explicitly to logical.
  obj$predicted_doublet <- toupper(as.character(scrub$predicted_doublet)) == "TRUE"

  obj$qc_flag <- if (!is.null(qc_flags[[srr]])) qc_flags[[srr]] else "none"

  n_pre <- ncol(obj)
  obj_filt <- subset(obj, subset = predicted_doublet == FALSE)
  n_post <- ncol(obj_filt)

  saveRDS(obj_filt, file.path("seurat_qc", paste0(srr, "_final.rds")))

  summary_rows[[srr]] <- data.frame(
    SRR = srr,
    n_cells_pre_doublet_filter = n_pre,
    n_doublets_removed = n_pre - n_post,
    n_cells_final = n_post,
    qc_flag = obj$qc_flag[1]
  )
  message(srr, ": ", n_pre - n_post, " doublets removed, ", n_post, " cells remain",
          if (obj$qc_flag[1] != "none") paste0(" [FLAGGED: ", obj$qc_flag[1], "]") else "")
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path("seurat_qc", "doublet_filter_summary.csv"), row.names = FALSE)
message("\n=== Stage 3 doublet filter summary ===")
print(summary_df)
message("\nStage 3 (Scrublet doublet removal) complete. Final per-sample objects: seurat_qc/<SRR>_final.rds")
