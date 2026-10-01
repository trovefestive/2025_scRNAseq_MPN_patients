#!/usr/bin/env Rscript
# Stage 7 (downstream-analysis-plan.md): cell type annotation.
#
# SingleR (celldex Human Primary Cell Atlas) at BOTH cluster and cell level, plus
# a curated MPN-BM marker panel (stage7_markers.R), at RNA_snn_res.0.8 (27 clusters).
#
# Also settles the four questions queued from Stages 4-6 (see notes):
#   Q1 cluster 6 -- real erythroid or junk? (holds ~76% of SRR27748964's cells)
#   Q2 cluster 9 -- 29.8% ribo + 584 nFeature: real low-RNA population or debris?
#   Q3 the SRR27748967-driven cluster -- drop it?
#   Q4 residual batch WITHIN each annotated cell type (Stage 5 left this open;
#      this is the test that separates leftover batch from real composition
#      difference, and gates any decision to tune Harmony theta)
#
# HPCA is read from the ExperimentHub cache under ~/.cache/R (populated on the
# login node) -- compute nodes have no outbound internet.
#
# Set STAGE7_TEST_N=<n> to dry-run on n cells per sample -> stage7_test/.

suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(BiocParallel); library(ggplot2); library(dplyr)
})

proj_dir <- Sys.getenv("PROJ_DIR", getwd())  # run from the project root
setwd(proj_dir)
options(future.globals.maxSize = 100 * 1024^3)
source("stage7_markers.R")

RES_COL <- "RNA_snn_res.0.8"
DIMS    <- 1:30
NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

test_n  <- Sys.getenv("STAGE7_TEST_N", "")
out_dir <- if (nzchar(test_n)) "stage7_test" else "stage7"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Loading stage6/clustered_umap.rds ===")
obj <- readRDS("stage6/clustered_umap.rds")

if (nzchar(test_n)) {
  n <- as.integer(test_n); set.seed(1)
  message("*** TEST MODE: ", n, " cells/sample -> ", out_dir, " ***")
  keep <- unlist(lapply(split(colnames(obj), obj$sample),
                        function(x) head(sample(x), min(n, length(x)))), use.names = FALSE)
  obj <- subset(obj, cells = keep)
}
Idents(obj) <- obj[[RES_COL]][, 1]
message("Cells: ", ncol(obj), " | clusters: ", nlevels(Idents(obj)))

# ---- SingleR --------------------------------------------------------------
message("\n=== Loading HPCA reference from cache ===")
ref <- HumanPrimaryCellAtlasData()
expr <- GetAssayData(obj, assay = "RNA", layer = "data")
bp <- MulticoreParam(workers = NCORES)

message("=== SingleR: cluster level (main + fine) ===")
sr_cl_main <- SingleR(test = expr, ref = ref, labels = ref$label.main,
                      clusters = Idents(obj), BPPARAM = bp)
sr_cl_fine <- SingleR(test = expr, ref = ref, labels = ref$label.fine,
                      clusters = Idents(obj), BPPARAM = bp)

message("=== SingleR: per cell (main) ===")
sr_cell <- SingleR(test = expr, ref = ref, labels = ref$label.main, BPPARAM = bp)
obj$singler_cell <- sr_cell$pruned.labels
obj$singler_cluster <- as.character(sr_cl_main$labels)[match(as.character(Idents(obj)),
                                                             rownames(sr_cl_main))]
saveRDS(list(cluster_main = sr_cl_main, cluster_fine = sr_cl_fine),
        file.path(out_dir, "singler_cluster_results.rds"))
message("SingleR done.")

# ---- Per-cluster consensus table ------------------------------------------
message("\n=== Building per-cluster annotation table ===")
sample_overall <- prop.table(table(obj$sample))
cl <- Idents(obj)

rows <- lapply(levels(cl), function(lv) {
  idx <- which(cl == lv)
  pc  <- table(obj$singler_cell[idx])
  pc  <- sort(pc[!is.na(names(pc))], decreasing = TRUE)
  sf  <- prop.table(table(obj$sample[idx]))
  data.frame(
    cluster            = lv,
    n_cells            = length(idx),
    pct_of_all         = round(100 * length(idx) / ncol(obj), 2),
    singler_cluster    = as.character(sr_cl_main$labels)[match(lv, rownames(sr_cl_main))],
    singler_fine       = as.character(sr_cl_fine$labels)[match(lv, rownames(sr_cl_fine))],
    singler_cell_top   = if (length(pc)) names(pc)[1] else NA_character_,
    singler_cell_purity= if (length(pc)) round(100 * pc[1] / sum(pc), 1) else NA_real_,
    median_nFeature    = round(median(obj$nFeature_RNA[idx])),
    median_pct_mt      = round(median(obj$percent.mt[idx]), 1),
    median_pct_ribo    = round(median(obj$percent.ribo[idx]), 1),
    top_sample         = names(which.max(sf)),
    max_enrichment     = round(max(sf / sample_overall[names(sf)]), 1),
    pct_from_SRR27748964 = round(100 * mean(obj$sample[idx] == "SRR27748964"), 2),
    pct_from_SRR27748967 = round(100 * mean(obj$sample[idx] == "SRR27748967"), 2),
    stringsAsFactors = FALSE)
})
ann <- do.call(rbind, rows)
ann <- ann[order(-ann$n_cells), ]
write.csv(ann, file.path(out_dir, "cluster_annotation.csv"), row.names = FALSE)
message("\n=== Per-cluster annotation ===")
print(ann[, c("cluster","n_cells","singler_cluster","singler_cell_top",
              "singler_cell_purity","median_nFeature","median_pct_ribo")], row.names = FALSE)

# ---- Marker expression per cluster ----------------------------------------
message("\n=== Curated marker mean expression per cluster ===")
present <- lapply(MARKERS, function(g) intersect(g, rownames(obj)))
avg <- AverageExpression(obj, features = unlist(present, use.names = FALSE),
                         assays = "RNA", layer = "data")$RNA
# AverageExpression prefixes numeric identities with "g" (g0, g1, ...). Strip it so
# cluster columns are indexable by the same labels used everywhere else in this script.
if (all(grepl("^g[0-9]+$", colnames(avg)))) colnames(avg) <- sub("^g", "", colnames(avg))
stopifnot(all(levels(cl) %in% colnames(avg)))
score <- t(sapply(names(present), function(k) {
  g <- present[[k]]
  if (!length(g)) return(rep(NA_real_, ncol(avg)))
  colMeans(log1p(avg[g, , drop = FALSE]))
}))
colnames(score) <- colnames(avg)
write.csv(round(score, 3), file.path(out_dir, "marker_scores_per_cluster.csv"))
top_sets <- apply(score, 2, function(x) paste(names(sort(x, decreasing = TRUE))[1:2], collapse = "/"))
message("\nTop 2 marker sets per cluster:")
print(data.frame(cluster = names(top_sets), top_marker_sets = unname(top_sets)), row.names = FALSE)

# ---- Q1-Q3: the specific clusters flagged at Stage 6 ----------------------
message("\n=== Q1/Q2/Q3: flagged clusters ===")
ery <- intersect(MARKERS$Erythroid, rownames(obj))
q_rows <- list()
for (lv in levels(cl)) {
  idx <- which(cl == lv)
  q_rows[[lv]] <- data.frame(
    cluster = lv, n_cells = length(idx),
    ery_score   = round(mean(log1p(avg[ery, lv])), 3),
    pct_HBB_pos = round(100 * mean(expr["HBB", idx] > 0), 1),
    median_nFeature = round(median(obj$nFeature_RNA[idx])),
    median_pct_mt   = round(median(obj$percent.mt[idx]), 1),
    median_pct_ribo = round(median(obj$percent.ribo[idx]), 1),
    singler = as.character(sr_cl_main$labels)[match(lv, rownames(sr_cl_main))])
}
qdf <- do.call(rbind, q_rows)
write.csv(qdf, file.path(out_dir, "erythroid_and_quality_check.csv"), row.names = FALSE)
message("\nErythroid score / quality by cluster (Q1 = cluster 6, Q2 = cluster 9):")
print(qdf[order(-qdf$ery_score), ], row.names = FALSE)

# ---- Q4: residual batch WITHIN each annotated cell type ------------------
# Stage 5 left dims 2-3 with high sample-linked variance and could not say whether
# that was leftover batch or genuine between-patient composition. Restricting to a
# single cell type removes composition as an explanation: whatever remains is batch.
message("\n=== Q4: batch variance within each SingleR cell type ===")
var_expl <- function(emb, g) apply(emb, 2, function(x) {
  gr <- sum((x - mean(x))^2); if (gr == 0) NA_real_ else 1 - sum((x - ave(x, g))^2) / gr })
emb_p <- Embeddings(obj, "pca")[, DIMS, drop = FALSE]
emb_h <- Embeddings(obj, "harmony")[, DIMS, drop = FALSE]
ct <- obj$singler_cell
keep_ct <- names(which(table(ct) >= 500))
b_rows <- lapply(keep_ct, function(k) {
  i <- which(ct == k); if (length(unique(obj$sample[i])) < 2) return(NULL)
  g <- droplevels(factor(obj$sample[i]))
  data.frame(cell_type = k, n_cells = length(i),
             n_samples = nlevels(g),
             r2_pca_pct     = round(100 * mean(var_expl(emb_p[i, , drop = FALSE], g), na.rm = TRUE), 2),
             r2_harmony_pct = round(100 * mean(var_expl(emb_h[i, , drop = FALSE], g), na.rm = TRUE), 2))
})
bdf <- do.call(rbind, b_rows)
bdf$reduction_pct <- round(100 * (bdf$r2_pca_pct - bdf$r2_harmony_pct) / bdf$r2_pca_pct, 1)
bdf <- bdf[order(-bdf$r2_harmony_pct), ]
write.csv(bdf, file.path(out_dir, "batch_within_celltype.csv"), row.names = FALSE)
print(bdf, row.names = FALSE)
message("\nInterpretation: high r2_harmony_pct WITHIN a cell type = residual batch ",
        "(composition cannot explain it). Low = Harmony succeeded and Stage 5's ",
        "dims 2-3 residual was real between-patient composition difference.")

# ---- Annotation-support markers ------------------------------------------
# Subsampled (500 cells/cluster) ON PURPOSE: this is to support annotation, it is
# NOT the Stage 9 DEG result, which runs on all cells with the plan's thresholds.
message("\n=== FindAllMarkers (max 500 cells/cluster -- annotation support only) ===")
mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,
                     max.cells.per.ident = 500, verbose = FALSE)
write.csv(mk, file.path(out_dir, "cluster_markers_subsampled.csv"), row.names = FALSE)
top10 <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")
write.csv(top10, file.path(out_dir, "cluster_top10_markers.csv"), row.names = FALSE)
message("\nTop 10 markers per cluster:")
print(as.data.frame(top10), row.names = FALSE)

# ---- Save + plots ---------------------------------------------------------
message("\n=== Saving ===")
saveRDS(obj, file.path(out_dir, "annotated.rds"))

ggsave(file.path(plot_dir, "umap_singler_cell.png"),
       DimPlot(obj, reduction = "umap", group.by = "singler_cell", label = TRUE,
               label.size = 3, repel = TRUE, raster = TRUE) + NoLegend() +
         ggtitle("SingleR (HPCA) per-cell label"), width = 9, height = 8, dpi = 150)
ggsave(file.path(plot_dir, "umap_singler_cluster.png"),
       DimPlot(obj, reduction = "umap", group.by = "singler_cluster", label = TRUE,
               label.size = 3, repel = TRUE, raster = TRUE) + NoLegend() +
         ggtitle("SingleR (HPCA) cluster-level label"), width = 9, height = 8, dpi = 150)
ggsave(file.path(plot_dir, "dotplot_curated_markers.png"),
       DotPlot(obj, features = present, cluster.idents = TRUE) +
         RotatedAxis() + theme(axis.text.x = element_text(size = 6)) +
         ggtitle("Curated MPN-BM markers by cluster"), width = 22, height = 8, dpi = 150)
ggsave(file.path(plot_dir, "umap_key_lineage_markers.png"),
       FeaturePlot(obj, reduction = "umap", raster = TRUE, ncol = 3,
                   features = c("CD34","HBB","PF4","LYZ","MS4A1","CD3D","NKG7","MPO","MKI67")),
       width = 13, height = 12, dpi = 150)

message("\n=== Stage 7 complete ===")
message("Object: ", file.path(out_dir, "annotated.rds"))
