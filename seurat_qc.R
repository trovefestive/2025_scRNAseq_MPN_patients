#!/usr/bin/env Rscript
# Stage 2 (downstream-analysis-plan.md): per-sample Seurat load + QC filtering.
# Filters: nFeature_RNA >= 200, percent.mt <= 20 (thresholds from plan, mirroring
# Jung et al. PMC11959246). Doublet removal (Scrublet) is Stage 3, done separately
# per sample BEFORE merge -- not handled here.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)

meta <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE) %>%
  select(Run, disease_stage, health_state, patient = Sample.Name, sex, age = AGE)

runs <- readLines("runs.txt")

out_dir <- "seurat_qc"
plot_dir <- file.path(out_dir, "qc_plots")
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

qc_rows <- list()

for (srr in runs) {
  message("=== Processing ", srr, " ===")
  mtx_dir <- file.path("cellranger_count", srr, "outs", "filtered_feature_bc_matrix")
  if (!dir.exists(mtx_dir)) {
    warning("Missing filtered_feature_bc_matrix/ for ", srr, " -- skipping")
    next
  }

  # Read10X (mtx dir) used instead of Read10X_h5 -- hdf5r is not installed in
  # scrnaseq-r and the mtx-format directory is equivalent, no h5 dependency needed.
  mat <- Read10X(mtx_dir)
  obj <- CreateSeuratObject(counts = mat, project = srr, min.cells = 3, min.features = 0)
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

  m <- meta %>% filter(Run == srr)
  stage_label <- if (nrow(m) == 0) NA_character_ else if (is.na(m$disease_stage) || m$disease_stage == "-") m$health_state else m$disease_stage

  obj$sample <- srr
  obj$patient <- if (nrow(m) > 0) m$patient else NA_character_
  obj$disease_stage <- stage_label
  obj$sex <- if (nrow(m) > 0) m$sex else NA_character_
  obj$age <- if (nrow(m) > 0) m$age else NA_character_

  n_pre <- ncol(obj)
  med_nfeature_pre <- median(obj$nFeature_RNA)
  med_ncount_pre <- median(obj$nCount_RNA)
  med_mt_pre <- median(obj$percent.mt)

  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) +
    plot_annotation(title = paste0(srr, " (", m$patient, ", ", stage_label, ") -- pre-filter, n=", n_pre))
  ggsave(file.path(plot_dir, paste0(srr, "_prefilter_violin.png")), p, width = 11, height = 4, dpi = 150)

  obj_filt <- subset(obj, subset = nFeature_RNA >= 200 & percent.mt <= 20)
  n_post <- ncol(obj_filt)

  saveRDS(obj_filt, file.path(out_dir, paste0(srr, ".rds")))

  qc_rows[[srr]] <- data.frame(
    SRR = srr,
    patient = m$patient,
    disease_stage = stage_label,
    n_cells_prefilter = n_pre,
    n_cells_postfilter = n_post,
    pct_retained = round(100 * n_post / n_pre, 1),
    median_nFeature_prefilter = med_nfeature_pre,
    median_nCount_prefilter = med_ncount_pre,
    median_percent_mt_prefilter = round(med_mt_pre, 2)
  )
}

summary_df <- bind_rows(qc_rows)
write.csv(summary_df, file.path(out_dir, "qc_summary.csv"), row.names = FALSE)

message("\n=== Stage 2 QC summary ===")
print(summary_df)
message("\nStage 2 (Seurat load + QC filter) complete. Per-sample filtered objects in ", out_dir, "/, summary at ", out_dir, "/qc_summary.csv")
