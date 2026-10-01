#!/usr/bin/env Rscript
# Stage 4 (downstream-analysis-plan.md): merge the 19 Stage-3 per-sample objects,
# LogNormalize + FindVariableFeatures(vst), ScaleData, PCA.
#
# Scope note: the plan lists Stage 4 as "Normalize" and puts PCA under Stage 6,
# but Harmony (Stage 5) takes a PCA reduction as input, so ScaleData + RunPCA run
# here to produce Stage 5's input. No integration and no clustering/UMAP happen in
# this script -- those stay Stages 5 and 6.
#
# Seurat v5 idiom: merge() keeps one counts layer per sample, so NormalizeData and
# FindVariableFeatures run per sample and HVGs are combined across samples (this is
# preferable to a global HVG call, which lets one deep sample dominate selection).
# Layers stay split for IntegrateLayers/Harmony in Stage 5.
#
# Set STAGE4_TEST_N=<n> to dry-run the whole code path on the first n samples,
# writing to stage4_test/ instead of stage4/.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)

options(future.globals.maxSize = 100 * 1024^3)

test_n <- Sys.getenv("STAGE4_TEST_N", "")
runs <- readLines("runs.txt")
if (nzchar(test_n)) {
  runs <- head(runs, as.integer(test_n))
  out_dir <- "stage4_test"
  message("*** TEST MODE: ", length(runs), " samples -> ", out_dir, " ***")
} else {
  out_dir <- "stage4"
}
dir.create(out_dir, showWarnings = FALSE)

# ---- Load Stage-3 objects -------------------------------------------------
message("=== Loading ", length(runs), " Stage-3 objects ===")
objs <- lapply(runs, function(srr) {
  p <- file.path("seurat_qc", paste0(srr, "_final.rds"))
  if (!file.exists(p)) stop("Missing Stage-3 object: ", p)
  o <- readRDS(p)
  message("  ", srr, ": ", ncol(o), " cells, ", nrow(o), " features")
  o
})
names(objs) <- runs

# ---- Merge ----------------------------------------------------------------
# add.cell.ids is required: 10x barcodes repeat across samples.
message("\n=== Merging ===")
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = runs, project = "MPN_BM")
rm(objs); invisible(gc())

message("Merged: ", nrow(merged), " features x ", ncol(merged), " cells")
message("Layers: ", length(Layers(merged[["RNA"]])), " (",
        paste(head(Layers(merged[["RNA"]]), 3), collapse = ", "), ", ...)")

stopifnot(ncol(merged) == sum(table(merged$sample)))
if (anyDuplicated(colnames(merged)) > 0) stop("Duplicate cell names after merge")

# ---- Normalize + HVG ------------------------------------------------------
message("\n=== NormalizeData (LogNormalize) ===")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 1e4,
                        verbose = FALSE)

message("=== FindVariableFeatures (vst, 2000) ===")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000,
                               verbose = FALSE)
hvg <- VariableFeatures(merged)
message("HVGs selected: ", length(hvg))
message("Top 20: ", paste(head(hvg, 20), collapse = ", "))

# ---- Scale + PCA ----------------------------------------------------------
message("\n=== ScaleData (HVGs only) ===")
merged <- ScaleData(merged, features = hvg, verbose = FALSE)

message("=== RunPCA (npcs=50) ===")
merged <- RunPCA(merged, features = hvg, npcs = 50, verbose = FALSE)

# ---- Save + diagnostics ---------------------------------------------------
message("\n=== Saving ===")
saveRDS(merged, file.path(out_dir, "merged_normalized.rds"))

write.csv(
  data.frame(
    sample = names(table(merged$sample)),
    n_cells = as.integer(table(merged$sample))
  ),
  file.path(out_dir, "stage4_cells_per_sample.csv"), row.names = FALSE
)
writeLines(hvg, file.path(out_dir, "variable_features.txt"))

ggsave(file.path(out_dir, "elbow_plot.png"),
       ElbowPlot(merged, ndims = 50) + ggtitle("Stage 4 PCA elbow (pre-integration)"),
       width = 7, height = 5, dpi = 150)

# Pre-integration PC1/PC2 by sample: the batch structure Harmony has to remove.
ggsave(file.path(out_dir, "pca_by_sample.png"),
       DimPlot(merged, reduction = "pca", group.by = "sample", raster = TRUE) +
         ggtitle("PC1/PC2 by sample -- PRE-integration"),
       width = 8, height = 6, dpi = 150)
ggsave(file.path(out_dir, "pca_by_stage.png"),
       DimPlot(merged, reduction = "pca", group.by = "disease_stage", raster = TRUE) +
         ggtitle("PC1/PC2 by disease stage -- PRE-integration"),
       width = 8, height = 6, dpi = 150)

sink(file.path(out_dir, "pca_loadings.txt"))
print(merged[["pca"]], dims = 1:10, nfeatures = 15)
sink()

message("\n=== Stage 4 complete ===")
message("Object: ", file.path(out_dir, "merged_normalized.rds"))
message("Cells: ", ncol(merged), " | Features: ", nrow(merged), " | HVGs: ", length(hvg))
print(table(merged$disease_stage))
