#!/usr/bin/env Rscript
# Stage 5 (downstream-analysis-plan.md): Harmony integration, batch = sample.
#
# Uses harmony::RunHarmony() directly rather than Seurat v5's IntegrateLayers()
# wrapper, because RunHarmony takes an explicit dims.use -- Stage 4's elbow put
# the knee at ~PC11-15, flat by 25-30, so we integrate over dims 1:30 rather than
# whatever npcs the wrapper defaults to. Layers are joined afterwards, which is
# what Stage 9's FindMarkers needs.
#
# Quantitative before/after: for each integrated dim, the fraction of variance
# explained by `sample` (one-way ANOVA R^2), on PCA vs Harmony embeddings. A big
# drop = Harmony removed sample-driven variance. This is the number to judge the
# stage on; the PC1/PC2 scatters are only the visual companion.
#
# Set STAGE5_TEST_N=<n> to dry-run on n cells per sample, writing to stage5_test/.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)

DIMS <- 1:30
BATCH <- "sample"

test_n <- Sys.getenv("STAGE5_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage5_test" else "stage5"
dir.create(out_dir, showWarnings = FALSE)

message("=== Loading stage4/merged_normalized.rds ===")
obj <- readRDS("stage4/merged_normalized.rds")
message("Loaded: ", nrow(obj), " features x ", ncol(obj), " cells")

if (nzchar(test_n)) {
  n <- as.integer(test_n)
  message("*** TEST MODE: subsampling ", n, " cells per sample -> ", out_dir, " ***")
  keep <- unlist(lapply(split(colnames(obj), obj$sample), function(cells) {
    head(sample(cells), min(n, length(cells)))
  }), use.names = FALSE)
  obj <- subset(obj, cells = keep)
  message("Subsampled to ", ncol(obj), " cells")
}

stopifnot("pca" %in% Reductions(obj))
stopifnot(BATCH %in% colnames(obj@meta.data))
message("Batch variable '", BATCH, "': ", length(unique(obj[[BATCH]][, 1])), " levels")

# ---- Harmony --------------------------------------------------------------
message("\n=== RunHarmony (dims ", min(DIMS), ":", max(DIMS), ", group.by = ", BATCH, ") ===")
obj <- RunHarmony(obj, group.by.vars = BATCH, reduction.use = "pca",
                  dims.use = DIMS, reduction.save = "harmony", verbose = TRUE)
message("harmony reduction: ", ncol(Embeddings(obj, "harmony")), " dims")

# ---- Batch-variance metric ------------------------------------------------
# R^2 of embedding ~ sample, per dim, before (pca) vs after (harmony).
var_expl <- function(emb, g) {
  apply(emb, 2, function(x) {
    grand <- sum((x - mean(x))^2)
    if (grand == 0) return(NA_real_)
    1 - sum((x - ave(x, g))^2) / grand
  })
}
g <- factor(obj[[BATCH]][, 1])
r2_pca <- var_expl(Embeddings(obj, "pca")[, DIMS, drop = FALSE], g)
r2_har <- var_expl(Embeddings(obj, "harmony")[, DIMS, drop = FALSE], g)

r2_df <- data.frame(dim = DIMS, r2_pca = round(r2_pca, 4), r2_harmony = round(r2_har, 4),
                    row.names = NULL)
r2_df$reduction_pct <- round(100 * (r2_df$r2_pca - r2_df$r2_harmony) / r2_df$r2_pca, 1)
write.csv(r2_df, file.path(out_dir, "batch_variance_explained.csv"), row.names = FALSE)

message("\n=== Variance explained by ", BATCH, " (mean over dims ",
        min(DIMS), ":", max(DIMS), ") ===")
message("  PCA     : ", round(100 * mean(r2_pca, na.rm = TRUE), 2), "%")
message("  Harmony : ", round(100 * mean(r2_har, na.rm = TRUE), 2), "%")
message("  Relative reduction: ",
        round(100 * (mean(r2_pca, na.rm = TRUE) - mean(r2_har, na.rm = TRUE)) /
                mean(r2_pca, na.rm = TRUE), 1), "%")
print(r2_df)

# ---- Join layers ----------------------------------------------------------
# Split-by-sample layers were needed for per-sample HVG selection and integration;
# downstream DE (Stage 9) needs a single counts/data layer.
message("\n=== JoinLayers ===")
message("Layers before: ", length(Layers(obj[["RNA"]])))
obj <- JoinLayers(obj)
message("Layers after : ", paste(Layers(obj[["RNA"]]), collapse = ", "))

# ---- Save + diagnostics ---------------------------------------------------
message("\n=== Saving ===")
saveRDS(obj, file.path(out_dir, "integrated_harmony.rds"))

ggsave(file.path(out_dir, "harmony_by_sample.png"),
       DimPlot(obj, reduction = "harmony", group.by = "sample", raster = TRUE) +
         ggtitle("Harmony_1/2 by sample -- POST-integration"),
       width = 8, height = 6, dpi = 150)
ggsave(file.path(out_dir, "harmony_by_stage.png"),
       DimPlot(obj, reduction = "harmony", group.by = "disease_stage", raster = TRUE) +
         ggtitle("Harmony_1/2 by disease stage -- POST-integration"),
       width = 8, height = 6, dpi = 150)

p <- ggplot(r2_df, aes(x = dim)) +
  geom_line(aes(y = 100 * r2_pca, colour = "PCA")) +
  geom_point(aes(y = 100 * r2_pca, colour = "PCA")) +
  geom_line(aes(y = 100 * r2_harmony, colour = "Harmony")) +
  geom_point(aes(y = 100 * r2_harmony, colour = "Harmony")) +
  labs(x = "dimension", y = paste0("% variance explained by ", BATCH),
       colour = NULL, title = "Batch variance per dimension, before vs after Harmony") +
  theme_classic()
ggsave(file.path(out_dir, "batch_variance_explained.png"), p, width = 8, height = 5, dpi = 150)

message("\n=== Stage 5 complete ===")
message("Object: ", file.path(out_dir, "integrated_harmony.rds"))
message("Cells: ", ncol(obj), " | Reductions: ", paste(Reductions(obj), collapse = ", "))
