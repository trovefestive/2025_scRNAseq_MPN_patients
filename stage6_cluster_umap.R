#!/usr/bin/env Rscript
# Stage 6 (downstream-analysis-plan.md): graph clustering + UMAP on the Harmony
# embedding from Stage 5.
#
# Clusters at a sweep of resolutions rather than one, and stores every result as a
# metadata column -- resolution is chosen at the Stage 6 check-in, not guessed here.
#
# Two QC questions carried in from earlier stages are answered per cluster:
#   1. Is any cluster ribosomal-driven? (Stage 4: PC6/7/10 were RPL*/RPS* axes)
#   2. Is any cluster driven by a single sample? Reported as an enrichment ratio
#      (cluster's fraction from a sample / that sample's fraction of all cells),
#      because raw fractions are misleading when samples range 878-13,930 cells.
#
# Set STAGE6_TEST_N=<n> to dry-run on n cells per sample, writing to stage6_test/.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)

DIMS <- 1:30
RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

test_n <- Sys.getenv("STAGE6_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage6_test" else "stage6"
dir.create(out_dir, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

message("=== Loading stage5/integrated_harmony.rds ===")
obj <- readRDS("stage5/integrated_harmony.rds")
message("Loaded: ", nrow(obj), " features x ", ncol(obj), " cells")

if (nzchar(test_n)) {
  n <- as.integer(test_n)
  message("*** TEST MODE: subsampling ", n, " cells per sample -> ", out_dir, " ***")
  set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$sample), function(cells) {
    head(sample(cells), min(n, length(cells)))
  }), use.names = FALSE)
  obj <- subset(obj, cells = keep)
  message("Subsampled to ", ncol(obj), " cells")
}

stopifnot("harmony" %in% Reductions(obj))

# Ribosomal fraction -- needed to answer QC question 1 above.
obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

# ---- Neighbours + clustering sweep ----------------------------------------
message("\n=== FindNeighbors (harmony, dims ", min(DIMS), ":", max(DIMS), ") ===")
obj <- FindNeighbors(obj, reduction = "harmony", dims = DIMS, verbose = FALSE)

message("=== FindClusters (resolutions: ", paste(RESOLUTIONS, collapse = ", "), ") ===")
obj <- FindClusters(obj, resolution = RESOLUTIONS, verbose = FALSE)

res_cols <- paste0("RNA_snn_res.", RESOLUTIONS)
n_clust <- sapply(res_cols, function(cc) nlevels(factor(obj@meta.data[[cc]])))
res_tbl <- data.frame(resolution = RESOLUTIONS, n_clusters = as.integer(n_clust),
                      row.names = NULL)
write.csv(res_tbl, file.path(out_dir, "resolution_sweep.csv"), row.names = FALSE)
message("\n=== Clusters per resolution ===")
print(res_tbl)

# ---- UMAP -----------------------------------------------------------------
message("\n=== RunUMAP (harmony, dims ", min(DIMS), ":", max(DIMS), ") ===")
obj <- RunUMAP(obj, reduction = "harmony", dims = DIMS,
               reduction.name = "umap", verbose = FALSE)

message("\n=== Saving object ===")
saveRDS(obj, file.path(out_dir, "clustered_umap.rds"))

# ---- Per-cluster QC across all resolutions --------------------------------
sample_overall <- prop.table(table(obj$sample))

qc_rows <- list()
for (i in seq_along(RESOLUTIONS)) {
  cc <- res_cols[i]
  cl <- factor(obj@meta.data[[cc]])
  for (lv in levels(cl)) {
    idx <- which(cl == lv)
    samp_frac <- prop.table(table(obj$sample[idx]))
    enrich <- samp_frac / sample_overall[names(samp_frac)]
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      resolution        = RESOLUTIONS[i],
      cluster           = lv,
      n_cells           = length(idx),
      pct_of_all_cells  = round(100 * length(idx) / ncol(obj), 2),
      n_samples_ge1pct  = sum(samp_frac >= 0.01),
      top_sample        = names(which.max(samp_frac)),
      top_sample_frac   = round(max(samp_frac), 3),
      max_enrichment    = round(max(enrich), 1),
      median_percent_ribo = round(median(obj$percent.ribo[idx]), 1),
      median_percent_mt   = round(median(obj$percent.mt[idx]), 1),
      median_nFeature     = round(median(obj$nFeature_RNA[idx])),
      pct_flagged_sample  = round(100 * mean(obj$qc_flag[idx] != "none"), 2)
    )
  }
}
qc_df <- do.call(rbind, qc_rows)
write.csv(qc_df, file.path(out_dir, "cluster_qc_all_resolutions.csv"), row.names = FALSE)

message("\n=== Potential problem clusters (any resolution) ===")
flagged <- qc_df[qc_df$median_percent_ribo >= 40 | qc_df$max_enrichment >= 5, ]
if (nrow(flagged) == 0) {
  message("None: no cluster has median percent.ribo >= 40% or a sample enriched >= 5x.")
} else {
  print(flagged)
}

# ---- Plots ----------------------------------------------------------------
message("\n=== Plotting ===")
for (i in seq_along(RESOLUTIONS)) {
  p <- DimPlot(obj, reduction = "umap", group.by = res_cols[i],
               label = TRUE, label.size = 4, raster = TRUE) +
    ggtitle(paste0("UMAP -- clusters at resolution ", RESOLUTIONS[i],
                   " (n=", res_tbl$n_clusters[i], ")")) + NoLegend()
  ggsave(file.path(plot_dir, paste0("umap_res", RESOLUTIONS[i], ".png")),
         p, width = 8, height = 7, dpi = 150)
}

for (v in c("sample", "disease_stage", "qc_flag")) {
  ggsave(file.path(plot_dir, paste0("umap_by_", v, ".png")),
         DimPlot(obj, reduction = "umap", group.by = v, raster = TRUE) +
           ggtitle(paste0("UMAP by ", v)),
         width = 9, height = 7, dpi = 150)
}

ggsave(file.path(plot_dir, "umap_qc_features.png"),
       FeaturePlot(obj, reduction = "umap",
                   features = c("percent.ribo", "percent.mt", "nFeature_RNA", "doublet_score"),
                   ncol = 2, raster = TRUE),
       width = 11, height = 9, dpi = 150)

message("\n=== Stage 6 complete ===")
message("Object: ", file.path(out_dir, "clustered_umap.rds"))
message("Cells: ", ncol(obj), " | Reductions: ", paste(Reductions(obj), collapse = ", "))
